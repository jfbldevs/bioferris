//! Shared create-only publication for trusted destination directories.

use std::fs::{self, OpenOptions};
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

#[derive(Debug)]
pub(crate) enum PublishError<E> {
    Io(io::Error),
    DestinationExists(PathBuf),
    TemporaryUnavailable,
    Write(E),
}

struct Temporary(PathBuf);

// Fault routing exists only in unit-test builds. Non-test expansions perform
// the original filesystem operation directly, without an injectable runtime API.
macro_rules! publication_io {
    ($point:ident, $operation:expr) => {{
        #[cfg(test)]
        {
            faults::perform(faults::Point::$point, || $operation)
        }
        #[cfg(not(test))]
        {
            $operation
        }
    }};
}

impl Drop for Temporary {
    fn drop(&mut self) {
        let _ = publication_io!(Cleanup, fs::remove_file(&self.0));
    }
}

// Publish only complete output using a same-directory hard link. No overwrite
// or copy fallback is allowed. This is not a hostile-directory or power-loss
// durability guarantee; callers must control the destination directory.
pub(crate) fn publish_new<E>(
    path: &Path,
    write: impl FnOnce(&mut dyn Write) -> Result<(), E>,
) -> Result<(), PublishError<E>> {
    static NEXT: AtomicU64 = AtomicU64::new(0);
    publish_with_counter(path, write, &NEXT)
}

fn publish_with_counter<E>(
    path: &Path,
    write: impl FnOnce(&mut dyn Write) -> Result<(), E>,
    counter: &AtomicU64,
) -> Result<(), PublishError<E>> {
    match publication_io!(Inspect, fs::symlink_metadata(path)) {
        Ok(_) => return Err(PublishError::DestinationExists(path.to_owned())),
        Err(error) if error.kind() == io::ErrorKind::NotFound => {}
        Err(error) => return Err(PublishError::Io(error)),
    }
    let parent = path
        .parent()
        .filter(|path| !path.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    let (temporary, mut file) = (0..128)
        .find_map(|_| {
            let candidate = parent.join(format!(
                ".seq-bio-{}-{}.tmp",
                std::process::id(),
                counter.fetch_add(1, Ordering::Relaxed)
            ));
            // The final destination must never double as our temporary file.
            // Both names use the same parent, including relative/symlink parents.
            if candidate.file_name() == path.file_name() {
                return None;
            }
            match publication_io!(
                Create,
                OpenOptions::new()
                    .write(true)
                    .create_new(true)
                    .open(&candidate)
            ) {
                Ok(file) => Some(Ok((Temporary(candidate), file))),
                Err(error) if error.kind() == io::ErrorKind::AlreadyExists => None,
                Err(error) => Some(Err(PublishError::Io(error))),
            }
        })
        .ok_or(PublishError::TemporaryUnavailable)??;
    #[cfg(test)]
    let result = faults::write_through(&mut file, write);
    #[cfg(not(test))]
    let result = write(&mut file);
    result.map_err(PublishError::Write)?;
    publication_io!(Sync, file.sync_all()).map_err(PublishError::Io)?;
    drop(file);
    publication_io!(Link, fs::hard_link(&temporary.0, path)).map_err(|error| {
        if error.kind() == io::ErrorKind::AlreadyExists {
            PublishError::DestinationExists(path.to_owned())
        } else {
            PublishError::Io(error)
        }
    })?;
    Ok(())
}

#[cfg(test)]
mod faults {
    use super::*;
    use std::cell::RefCell;

    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    pub(super) enum Point {
        Inspect,
        Create,
        Flush,
        Sync,
        Link,
        Cleanup,
    }

    #[derive(Default)]
    pub(super) struct Plan {
        pub(super) failure: Option<(Point, usize, io::ErrorKind)>,
        pub(super) write_after: Option<usize>,
        pub(super) max_write: Option<usize>,
        pub(super) interrupted: bool,
        pub(super) zero_write: bool,
        pub(super) calls: [usize; 6],
        pub(super) written: usize,
    }

    thread_local! {
        static ACTIVE: RefCell<Option<Plan>> = const { RefCell::new(None) };
    }

    pub(super) struct Guard {
        owner: std::thread::ThreadId,
    }

    impl Guard {
        pub(super) fn install(plan: Plan) -> Self {
            ACTIVE.with(|active| {
                assert!(active.borrow().is_none(), "fault plans cannot be nested");
                *active.borrow_mut() = Some(plan);
            });
            Self {
                owner: std::thread::current().id(),
            }
        }

        pub(super) fn calls(&self, point: Point) -> usize {
            assert_eq!(self.owner, std::thread::current().id());
            ACTIVE.with(|active| active.borrow().as_ref().unwrap().calls[point as usize])
        }
    }

    impl Drop for Guard {
        fn drop(&mut self) {
            assert_eq!(self.owner, std::thread::current().id());
            ACTIVE.with(|active| *active.borrow_mut() = None);
        }
    }

    pub(super) fn perform<T>(
        point: Point,
        operation: impl FnOnce() -> io::Result<T>,
    ) -> io::Result<T> {
        let failure = ACTIVE.with(|active| {
            let mut active = active.borrow_mut();
            let plan = active.as_mut()?;
            plan.calls[point as usize] += 1;
            plan.failure.and_then(|(target, occurrence, kind)| {
                (point == target && plan.calls[point as usize] == occurrence).then_some(kind)
            })
        });
        if let Some(kind) = failure {
            Err(io::Error::new(kind, "injected publication fault"))
        } else {
            operation()
        }
    }

    struct Writer<'a>(&'a mut dyn Write);

    impl Write for Writer<'_> {
        fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
            let limit = ACTIVE.with(|active| {
                let mut active = active.borrow_mut();
                let Some(plan) = active.as_mut() else {
                    return Ok(bytes.len());
                };
                if std::mem::take(&mut plan.interrupted) {
                    return Err(io::Error::new(
                        io::ErrorKind::Interrupted,
                        "retry this write",
                    ));
                }
                if plan.zero_write {
                    return Ok(0);
                }
                let remaining = plan
                    .write_after
                    .map_or(usize::MAX, |after| after.saturating_sub(plan.written));
                if remaining == 0 && !bytes.is_empty() {
                    return Err(io::Error::other("injected publication write failure"));
                }
                Ok(bytes
                    .len()
                    .min(remaining)
                    .min(plan.max_write.unwrap_or(usize::MAX)))
            })?;
            let written = self.0.write(&bytes[..limit])?;
            ACTIVE.with(|active| {
                if let Some(plan) = active.borrow_mut().as_mut() {
                    plan.written += written;
                }
            });
            Ok(written)
        }

        fn flush(&mut self) -> io::Result<()> {
            perform(Point::Flush, || self.0.flush())
        }
    }

    pub(super) fn write_through<E>(
        writer: &mut dyn Write,
        callback: impl FnOnce(&mut dyn Write) -> Result<(), E>,
    ) -> Result<(), E> {
        callback(&mut Writer(writer))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    struct Directory(PathBuf);

    impl Directory {
        fn new() -> Self {
            static NEXT: AtomicU64 = AtomicU64::new(0);
            let path = std::env::temp_dir().join(format!(
                "bioferris-publication-test-{}-{}",
                std::process::id(),
                NEXT.fetch_add(1, Ordering::Relaxed)
            ));
            fs::create_dir(&path).unwrap();
            Self(path)
        }
    }

    impl Drop for Directory {
        fn drop(&mut self) {
            fs::remove_dir_all(&self.0).unwrap();
        }
    }

    #[test]
    fn destination_matching_temporary_name_is_never_written_before_publication() {
        let directory = Directory::new();
        let destination = directory
            .0
            .join(format!(".seq-bio-{}-0.tmp", std::process::id()));
        publish_with_counter(
            &destination,
            |writer| -> io::Result<()> {
                assert!(!destination.exists());
                writer.write_all(b"complete result")?;
                assert!(!destination.exists());
                Ok(())
            },
            &AtomicU64::new(0),
        )
        .unwrap();
        assert_eq!(fs::read(&destination).unwrap(), b"complete result");
        assert_eq!(fs::read_dir(&directory.0).unwrap().count(), 1);
    }

    #[test]
    fn failed_callback_removes_temporary_data_without_publishing() {
        let directory = Directory::new();
        let destination = directory.0.join("result.txt");
        let result = publish_new(&destination, |writer| -> io::Result<()> {
            writer.write_all(b"incomplete")?;
            Err(io::Error::other("deliberate callback failure"))
        });
        assert!(matches!(result, Err(PublishError::Write(_))));
        assert!(!destination.exists());
        assert_eq!(fs::read_dir(&directory.0).unwrap().count(), 0);
    }

    #[test]
    fn destination_created_during_write_is_preserved_and_reported() {
        let directory = Directory::new();
        let destination = directory.0.join("result.txt");
        let result = publish_new(&destination, |writer| -> io::Result<()> {
            writer.write_all(b"our result")?;
            fs::write(&destination, b"other result")?;
            Ok(())
        });
        assert!(
            matches!(result, Err(PublishError::DestinationExists(path)) if path == destination)
        );
        assert_eq!(fs::read(&destination).unwrap(), b"other result");
        assert_eq!(fs::read_dir(&directory.0).unwrap().count(), 1);
    }

    #[test]
    fn filesystem_failures_stop_at_the_real_publication_boundary() {
        use faults::{Guard, Plan, Point};
        for point in [Point::Inspect, Point::Create, Point::Sync, Point::Link] {
            let directory = Directory::new();
            let destination = directory.0.join("result.txt");
            let fault = Guard::install(Plan {
                failure: Some((point, 1, io::ErrorKind::PermissionDenied)),
                ..Plan::default()
            });
            let mut called = false;
            let result = publish_new(&destination, |writer| -> io::Result<()> {
                called = true;
                writer.write_all(b"complete result")
            });
            assert!(matches!(result, Err(PublishError::Io(error))
                if error.kind() == io::ErrorKind::PermissionDenied
                    && error.to_string() == "injected publication fault"));
            assert_eq!(called, matches!(point, Point::Sync | Point::Link));
            assert_eq!(fault.calls(point), 1);
            assert_eq!(fault.calls(Point::Cleanup), usize::from(called));
            assert!(!destination.exists());
            assert_eq!(fs::read_dir(&directory.0).unwrap().count(), 0);
        }
    }

    #[test]
    fn every_partial_write_failure_removes_real_temporary_data() {
        use faults::{Guard, Plan, Point};
        const PAYLOAD: &[u8] = b">dna\nATGCatgcNN\n";
        for prefix in 0..PAYLOAD.len() {
            let directory = Directory::new();
            let destination = directory.0.join("result.fa");
            let fault = Guard::install(Plan {
                write_after: Some(prefix),
                max_write: Some(3),
                ..Plan::default()
            });
            let result = publish_new(&destination, |writer| writer.write_all(PAYLOAD));
            assert!(matches!(result, Err(PublishError::Write(error))
                if error.to_string() == "injected publication write failure"));
            assert_eq!(fault.calls(Point::Sync), 0);
            assert_eq!(fault.calls(Point::Link), 0);
            assert_eq!(fault.calls(Point::Cleanup), 1);
            assert!(!destination.exists());
            assert_eq!(fs::read_dir(&directory.0).unwrap().count(), 0);
        }
    }

    #[test]
    fn interrupted_and_short_writes_preserve_the_complete_payload() {
        let directory = Directory::new();
        let destination = directory.0.join("result.txt");
        let fault = faults::Guard::install(faults::Plan {
            interrupted: true,
            max_write: Some(1),
            ..faults::Plan::default()
        });
        publish_new(&destination, |writer| -> io::Result<()> {
            writer.write_all(b"every byte survives")?;
            writer.flush()
        })
        .unwrap();
        assert_eq!(fs::read(&destination).unwrap(), b"every byte survives");
        assert_eq!(fault.calls(faults::Point::Flush), 1);
        assert_eq!(fault.calls(faults::Point::Sync), 1);
        assert_eq!(fault.calls(faults::Point::Link), 1);
        assert_eq!(fault.calls(faults::Point::Cleanup), 1);
        assert_eq!(fs::read_dir(&directory.0).unwrap().count(), 1);
    }

    #[test]
    fn zero_write_and_explicit_callback_flush_failure_do_not_publish() {
        for flush in [false, true] {
            let directory = Directory::new();
            let destination = directory.0.join("result.txt");
            let _fault = faults::Guard::install(faults::Plan {
                failure: flush.then_some((faults::Point::Flush, 1, io::ErrorKind::Other)),
                zero_write: !flush,
                ..faults::Plan::default()
            });
            let result = publish_new(&destination, |writer| -> io::Result<()> {
                writer.write_all(b"result")?;
                writer.flush()
            });
            assert!(matches!(result, Err(PublishError::Write(error))
                if error.kind() == if flush { io::ErrorKind::Other } else { io::ErrorKind::WriteZero }));
            assert!(!destination.exists());
            assert_eq!(fs::read_dir(&directory.0).unwrap().count(), 0);
        }
    }

    #[test]
    fn callback_error_type_and_panic_are_not_replaced_by_publication() {
        #[derive(Debug, PartialEq)]
        struct CallbackError(u32);

        let directory = Directory::new();
        let destination = directory.0.join("result.txt");
        let result = publish_new(&destination, |writer| {
            writer.write_all(b"private").unwrap();
            Err(CallbackError(42))
        });
        assert!(matches!(
            result,
            Err(PublishError::Write(CallbackError(42)))
        ));
        let panic = std::panic::catch_unwind(|| {
            let _ = publish_new(&destination, |writer| -> io::Result<()> {
                writer.write_all(b"private")?;
                panic!("callback panic");
            });
        });
        assert!(panic.is_err());
        assert!(!destination.exists());
        assert_eq!(fs::read_dir(&directory.0).unwrap().count(), 0);
    }

    #[test]
    fn cleanup_failure_is_best_effort_and_preserves_the_original_result() {
        for callback_succeeds in [false, true] {
            let directory = Directory::new();
            let destination = directory.0.join("result.txt");
            let fault = faults::Guard::install(faults::Plan {
                failure: Some((faults::Point::Cleanup, 1, io::ErrorKind::PermissionDenied)),
                ..faults::Plan::default()
            });
            let result = publish_new(&destination, |writer| -> io::Result<()> {
                writer.write_all(b"private or complete")?;
                if callback_succeeds {
                    Ok(())
                } else {
                    Err(io::Error::other("original callback error"))
                }
            });
            if callback_succeeds {
                result.unwrap();
                assert_eq!(fs::read(&destination).unwrap(), b"private or complete");
            } else {
                assert!(matches!(result, Err(PublishError::Write(error))
                    if error.to_string() == "original callback error"));
                assert!(!destination.exists());
            }
            assert_eq!(fault.calls(faults::Point::Cleanup), 1);
            // This is the existing best-effort Drop contract, not a guarantee of
            // removal when the filesystem refuses it. Directory owns all files.
            assert_eq!(
                fs::read_dir(&directory.0).unwrap().count(),
                if callback_succeeds { 2 } else { 1 }
            );
        }
    }

    #[test]
    fn injected_link_collision_uses_destination_exists_classification() {
        let directory = Directory::new();
        let destination = directory.0.join("result.txt");
        let _fault = faults::Guard::install(faults::Plan {
            failure: Some((faults::Point::Link, 1, io::ErrorKind::AlreadyExists)),
            ..faults::Plan::default()
        });
        let result = publish_new(&destination, |writer| writer.write_all(b"result"));
        assert!(
            matches!(result, Err(PublishError::DestinationExists(path)) if path == destination)
        );
        assert_eq!(fs::read_dir(&directory.0).unwrap().count(), 0);
    }

    #[test]
    fn exhausting_all_temporary_candidates_never_truncates_existing_files() {
        let directory = Directory::new();
        let destination = directory.0.join("result.txt");
        for number in 0..128 {
            fs::write(
                directory
                    .0
                    .join(format!(".seq-bio-{}-{number}.tmp", std::process::id())),
                format!("sentinel {number}"),
            )
            .unwrap();
        }
        let result = publish_with_counter(
            &destination,
            |_| -> io::Result<()> { panic!("callback must not run") },
            &AtomicU64::new(0),
        );
        assert!(matches!(result, Err(PublishError::TemporaryUnavailable)));
        assert!(!destination.exists());
        assert_eq!(fs::read_dir(&directory.0).unwrap().count(), 128);
        for number in 0..128 {
            assert_eq!(
                fs::read_to_string(
                    directory
                        .0
                        .join(format!(".seq-bio-{}-{number}.tmp", std::process::id()))
                )
                .unwrap(),
                format!("sentinel {number}")
            );
        }
    }

    #[test]
    fn occupied_file_and_directory_are_rejected_before_callback() {
        let directory = Directory::new();
        for is_directory in [false, true] {
            let destination = directory
                .0
                .join(if is_directory { "directory" } else { "file" });
            if is_directory {
                fs::create_dir(&destination).unwrap();
            } else {
                fs::write(&destination, b"sentinel").unwrap();
            }
            let result = publish_new(&destination, |_| -> io::Result<()> {
                panic!("occupied destination must not invoke callback")
            });
            assert!(
                matches!(result, Err(PublishError::DestinationExists(path)) if path == destination)
            );
            if !is_directory {
                assert_eq!(fs::read(&destination).unwrap(), b"sentinel");
            }
        }
        assert_eq!(fs::read_dir(&directory.0).unwrap().count(), 2);
    }

    #[cfg(unix)]
    #[test]
    fn existing_and_dangling_symlinks_are_not_followed_or_replaced() {
        use std::os::unix::fs::symlink;
        let directory = Directory::new();
        let target = directory.0.join("target");
        fs::write(&target, b"sentinel").unwrap();
        for name in ["target", "missing"] {
            let destination = directory.0.join(format!("link-to-{name}"));
            let expected_target = directory.0.join(name);
            symlink(&expected_target, &destination).unwrap();
            let result = publish_new(&destination, |_| -> io::Result<()> {
                panic!("symlink destination must not invoke callback")
            });
            assert!(
                matches!(result, Err(PublishError::DestinationExists(path)) if path == destination)
            );
            assert_eq!(fs::read_link(&destination).unwrap(), expected_target);
        }
        assert_eq!(fs::read(&target).unwrap(), b"sentinel");
        assert!(!directory.0.join("missing").exists());
    }

    #[test]
    fn competing_publishers_have_exactly_one_complete_winner() {
        use std::sync::Arc;
        use std::sync::atomic::AtomicUsize;
        use std::time::{Duration, Instant};

        let directory = Directory::new();
        let destination = directory.0.join("winner.txt");
        let arrived = Arc::new(AtomicUsize::new(0));
        let results = std::thread::scope(|scope| {
            let workers = (0..4)
                .map(|worker| {
                    let arrived = Arc::clone(&arrived);
                    let destination = &destination;
                    scope.spawn(move || {
                        let payload = format!("complete payload from publisher {worker}\n");
                        let result = publish_new(destination, |writer| -> io::Result<()> {
                            writer.write_all(payload.as_bytes())?;
                            arrived.fetch_add(1, Ordering::SeqCst);
                            let deadline = Instant::now() + Duration::from_secs(5);
                            while arrived.load(Ordering::SeqCst) != 4 {
                                if Instant::now() >= deadline {
                                    return Err(io::Error::new(
                                        io::ErrorKind::TimedOut,
                                        "publisher rendezvous deadline",
                                    ));
                                }
                                std::thread::yield_now();
                            }
                            Ok(())
                        });
                        (payload, result)
                    })
                })
                .collect::<Vec<_>>();
            workers
                .into_iter()
                .map(|worker| worker.join().unwrap())
                .collect::<Vec<_>>()
        });
        assert_eq!(
            results.iter().filter(|(_, result)| result.is_ok()).count(),
            1
        );
        for (_, result) in &results {
            assert!(
                result.is_ok()
                    || matches!(result,
                Err(PublishError::DestinationExists(path)) if path == &destination)
            );
        }
        let winner = results.iter().find(|(_, result)| result.is_ok()).unwrap();
        assert_eq!(fs::read(&destination).unwrap(), winner.0.as_bytes());
        assert_eq!(fs::read_dir(&directory.0).unwrap().count(), 1);
    }
}
