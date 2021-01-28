# Safe-Error-Attacks-on-SIKE-and-CSIDH

This repository contains auxiliary material for the paper "Safe-Error Attacks on SIKE and CSIDH".

Authors
 - [Fabio Campos](https://www.sopmac.de/) `<campos@sopmac.de>` 
 - Marcel Müller `<marcel@qpc.tu-darmstadt.de>`
 - Juliane Krämmer `<juliane@qpc.tu-darmstadt.de>`

<!-- The paper is available at https://eprint.iacr.org/2020/1005.pdf -->

# Overview

This archive contains the following 
- `CSIDH` contains a slightly modified CSIDH implementation from [here](https://github.com/csidhfi/csidhfi)

- `SIKEp434` contains a slightly modified CSIDH implementation from [here](https://github.com/mupq/pqm4)
  
- `attack` contains the attack scripts used in Section 5 of the paper. If you have a ChipWhisperer-Lite (CW1173) 32-bit basic board which features a STM32F303, feel free to reproduce the results presented in the paper.


# Licenses

Code in this repository that does not indicate otherwise is placed in the public domain. 

For the third party code see their licenses:
- [pqm4](https://github.com/mupq/pqm4): [https://creativecommons.org/publicdomain/zero/1.0/](https://creativecommons.org/publicdomain/zero/1.0/)
- [libopencm3](https://github.com/libopencm3/libopencm3): [https://github.com/libopencm3/libopencm3/blob/master/COPYING.LGPL3](https://github.com/newaetech/chipwhisperer/blob/develop/LICENSE.txt)
- [chipwhisperer](https://github.com/newaetech/chipwhisperer): [https://github.com/newaetech/chipwhisperer/blob/develop/LICENSE.txt](https://github.com/newaetech/chipwhisperer/blob/develop/LICENSE.txt)