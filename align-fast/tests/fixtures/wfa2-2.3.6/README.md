# WFA2 2.3.6 affine fixture

`wavefront.json` contains 19 immutable exact global gap-affine observations generated
by the isolated reference-validation laboratory. It invokes the documented
`align_benchmark` CLI from official WFA2-lib `v2.3.6` at commit
`bcf473a6561fa297934a80eaa7c04b4ee412360c`, using high-memory traceback, one thread,
and no heuristic.

The external generator independently replays every reported CIGAR and checks its
penalty. Only this inert JSON enters the Rust repository; WFA2 source, binaries,
headers, and build tools remain outside the workspace and are never linked,
imported, or invoked by Cargo.

The reference implementation is MIT licensed and is used solely as an experimental
oracle and performance comparator.
