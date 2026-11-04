# Third-party notices

## Exponential Ziggurat sampler

The exponential tables and algorithm in `crates/crn-craft/src/ssa/ziggurat*`
are adapted from Julia 1.13.3, `stdlib/Random/src/normal.jl`. The adaptation
uses CRN-craft random streams, explicit work/cancellation checks, iterative
rejection and strict positive-wait checks; it does not promise Julia RNG replay.

Upstream: https://github.com/JuliaLang/julia/blob/v1.13.3/stdlib/Random/src/normal.jl

MIT License

Copyright (c) 2009-2025: Jeff Bezanson, Stefan Karpinski, Viral B. Shah, and other contributors: https://github.com/JuliaLang/julia/contributors

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

end of terms and conditions

[Julia's LICENSE.md ends with a pointer to its own THIRDPARTY.md, which
covers other software in the Julia repository; none of it is used here.]

## Optional parameter-optimization references

The isolated `validation/parameter-optimization/worker` workspace uses rebop
0.10.0 as an optional benchmark reference. Its generated wrappers invoke the
upstream macro API; it is not a dependency or implementation of the default
CRN-craft engine. The upstream MIT notice is preserved in
`validation/parameter-optimization/REBOP_LICENSE.txt`. Cargo.lock pins its
reference dependencies. Python fcmaes-rust/rebop dependencies are separately
listed in that directory's requirements file.
