# MAFFT 7.526 black-box observations

`fft-ns-1-default-v1.json`, `fft-ns-2-default-v1.json`,
`fft-ns-i-default-v1.json`, `l-ins-i-localpair-v1.json`,
`l-ins-i-offset-sign-v1.json`, `e-ins-i-oldgenafpair-v1.json`, and
`e-ins-i-parameter-sweep-v1.json` are inert
observation records produced outside the Rust workspace by the official Linux
amd64 package, invoked respectively as:

```text
mafft --retree 1 --maxiterate 0 --thread 1 --quiet INPUT
mafft --retree 2 --maxiterate 0 --thread 1 --inputorder --quiet INPUT
mafft --retree 2 --maxiterate 1000 --thread 1 --inputorder --quiet INPUT
mafft --localpair --maxiterate 1000 --thread 1 --inputorder --quiet INPUT
mafft --localpair --maxiterate 0 --thread 1 --inputorder --treeout --quiet
      --lep VALUE INPUT
mafft --genafpair --oldgenafpair --ep 0 --maxiterate 1000 --thread 1
      --inputorder --quiet INPUT
mafft --genafpair --oldgenafpair --ep 0 --maxiterate 0 --thread 1
      --inputorder --treeout --quiet [OPTION VALUE] INPUT
```

`g-ins-i-globalpair-v1.json` was produced in an isolated temporary laboratory
with the official signed Windows package and this opaque invocation:

```text
mafft.bat --globalpair --maxiterate 1000 --thread 1 --inputorder --quiet
          --out OUTPUT INPUT
```

Package SHA-256:
`b8298e379bcd244be0e12f325461d863b93571da78c4499c2a2e4c252679c9de`.
G-INS-i fixture SHA-256:
`7d0e23625119606a592e25b3113eb1d92a15caba180c0b6b9963eb14fb846906`.
L-INS-i fixture SHA-256:
`4c1273b81a396c19670094ba03c3a63ca3c8613b042345d435cfd838d45ea55a`.
E-INS-i alignment fixture SHA-256:
`20d0be6db389901b27c9482bfe7ba7aa73d75ac51556375dc103e7172e2e0af8`.
The compact E-INS-i parameter fixture SHA-256 is
`f06fc949998ea77dbd2c91ebcaa46d862f2ac91c53d64f4e06194ae05da40991`;
it summarizes the isolated raw observation with SHA-256
`93bc2c55de6bf1a418f4adbf6308b70850b065bf320206f024113c79ff3fc700`.
The compact L-INS-i offset-sign fixture summarizes an isolated raw observation
record with SHA-256
`9bccaa57c0e623f929f1675fed4520020ea3f22af79ba8952eac473182884f8c`.
FFT-NS-2 fixture SHA-256:
`96805cb0a35a1e22a73b14242507b8fd6323aba5f84f23c5136eef1ff91dfd25`.
FFT-NS-i fixture SHA-256:
`4c9900c326eb5f4c0fa0d3dc52304d18a53bfc57894a4cb8a20c421d2757ce78`.
Official Windows package SHA-256:
`7c3a508c8512af04007d68019c6fbd4ff00d891244b58f87ec9aa0e53fbae8a7`.

The generator, input specification, official package, and rerun instructions
live in the sibling isolated validation laboratory. No MAFFT source was
downloaded or inspected. This library does not invoke the executable.

The Rust FFT profiles target the original 2002 publication, G/L-INS-i target
ancestral MAFFT from the 2005 publications, and E-INS-i targets its explicit
2008 description plus Altschul's 1998 generalized affine model. MAFFT 7.526
has later defaults and heuristics, recorded in the fixtures, so nontrivial
comparisons use SP and TC rather than requiring byte-identical output.
