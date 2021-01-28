/*
    This file is part of the ChipWhisperer Example Targets
    Copyright (C) 2012-2017 NewAE Technology Inc.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "hal.h"
#include "simpleserial.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "fp.h"
#include "csidh.h"
#include "mont.h"
#include "uint.h"

#define REDUCEDSK

uint8_t csidh_glitch(uint8_t *in)
{

//#define KEY_0
    bool error = 0;
#ifdef KEY_0

    // private_key priv_alice = {{0, 1}};
    public_key shared_secret = {{{0xc3c2a07e486a9124, 0x162c9e08091fcdae, 0x6256ba93cc79389a, 0xf9d869529eae7626, 0xe308cf9a9e61fd76, 0xc5573bb4df3b9f6d, 0xb0ca1c2e87878bba, 0x1c73e3ab70aa49a2}}};

    proj A = {{{0xa7071cf2062c5b28, 0x4ef6c4e374631ad5, 0x075a4dd6d3013833, 0xa3c0a67a26b9e943, 0x51601a8a437952f2, 0x4f45902681f6b516, 0xc8364e54abb888f4, 0x0266ebb102f36783}}, {{0xc8fc8df598726f0a, 0x7b1bc81750a6af95, 0x5d319e67c1e961b4, 0xb0aa7275301955f1, 0x4a080672d9ba6c64, 0x97a5ef8a246ee77b, 0x06ea9e5d4383676a, 0x3496e2e117e0ec80}}};
    proj P = {{{0x3a035b56ffcf85a1, 0x26c52e62d79cd414, 0xdf35a4f0b45f71cb, 0xd291c8ff2ef03951, 0x825ec76038b3abca, 0xcc79e0c1365f9d7f, 0x87330014e1ecacb8, 0x1b5d3aebb93e3c89}}, {{0xc5aa1148b5605995, 0xd692a34936d5e7ec, 0x0f51fad0c047b6ff, 0x091aa823ebb7480a, 0xe6e79d572eb71c92, 0xf1c1ea7a99b44398, 0x5de7bf363594f415, 0x1c914b27d040e414}}};
    proj Pd = {{{0xc96f63361f4096b5, 0x3e53ebf8133e0279, 0x38c84a1b06cd0c91, 0x8de55edaebc6fc48, 0x4d8cb007d9fb7e6f, 0xd77c82ab93436d27, 0x7f3aa7e4aaa30313, 0x418e614a17616a19}}, {{0x3411a57153b97222, 0xe0e8ad126a875bed, 0x968bc02a1b267f9d, 0x6a1a715258f1a203, 0x9230180c312c2f23, 0x47104975bc2dba23, 0x30eb8a0edeb0b4c3, 0x55e5d37a71904e1f}}};
    proj K = {{{0xe078674610af7eaf, 0x92f8be8e112d8644, 0xdccf58317f042d0b, 0x324c9301302ede9e, 0x00275a83552d5f6c, 0xa460771fcc5a9f28, 0x5fab6ae359e5ac83, 0x288d6d70ecdafbbb}}, {{0xf7bed5a5950e83ad, 0x8099f93da24c3802, 0x7c07d9c5feabbeec, 0x71905cff8bd76fcb, 0x17977840811116ad, 0x566de9ad4ac853c9, 0x57ae3db851d07da3, 0x4ab797632d0a1ac5}}};

    // dummy
    error |= xISOG(&A, &P, &Pd, &K, 359, 1);

    fp_cswap(&P.x, &Pd.x, 0);
    fp_cswap(&P.z, &Pd.z, 0);
    uint_c cof = {{1}};
    uint_c l;
    xMUL(&K, &A, &P, &cof);
    uint_set(&l, 353);
    xMUL(&Pd, &A, &Pd, &l);

    // real
    error |= xISOG(&A, &P, &Pd, &K, 353, 0);

    fp_inv(&A.z);
    fp_mul2(&A.x, &A.z);

#else

    // private_key priv_alice = {{-1, 1}};
    public_key shared_secret = {{{0x40e3b5008005eb86, 0x49b705e6d9b47265, 0xa6786ef78224d926, 0x4323e4d98dd2f44d, 0x0a651bbc99f0a1de, 0x32df18c48eb96a72, 0x36a5dcda44dd9b65, 0x150cd06495afd0f5}}};

    proj A = {{{0xa7071cf2062c5b28, 0x4ef6c4e374631ad5, 0x075a4dd6d3013833, 0xa3c0a67a26b9e943, 0x51601a8a437952f2, 0x4f45902681f6b516, 0xc8364e54abb888f4, 0x0266ebb102f36783}}, {{0xc8fc8df598726f0a, 0x7b1bc81750a6af95, 0x5d319e67c1e961b4, 0xb0aa7275301955f1, 0x4a080672d9ba6c64, 0x97a5ef8a246ee77b, 0x06ea9e5d4383676a, 0x3496e2e117e0ec80}}};
    proj P = {{{0xcec97b1b06411f36, 0x8f1670ebc825a151, 0xb320c7a44560f223, 0xb5caa5e261bed15d, 0x035668fd252c2d14, 0x72e043a0f7299705, 0x3b758e031fae7245, 0x54e21834458e8e79}}, {{0xa98414332db5e44f, 0x3cde5baa45a53421, 0x8579f3bf2895b4f4, 0x5d6a3f90e096e9db, 0xd66b9ec9eb44a080, 0x2e22d1f41c929e52, 0x4bf4f9eba833b7c8, 0x054a72f57ad1e8b7}}};
    proj Pd = {{{0x34cc626f74594a9f, 0x98d571668026e7cf, 0x2cec58fb51704b85, 0x1b1ef60cbdf0342f, 0x1fc7bf32f1ac0021, 0x3372e6e76c855fe3, 0x01e73078eb21507e, 0x3e65e17421a361b9}}, {{0x0200d33a55afb57a, 0xe7543c8b412a44f5, 0x43c78d2cf2bdb1fd, 0x60ac555d8e2b3f5a, 0x1fc7cf0e78500ed3, 0x28d285f83a060858, 0xfa76e8a7ec2b49c0, 0x060d692e52acb29f}}};
    proj K = {{{0x85aa6c62d1068ea5, 0x1b21fddf7d415553, 0xdc95c09bf1b1af38, 0x132325aa3ac3626f, 0xcac04212a2d2e974, 0x3fcdc1fa8988be5f, 0x0177232dd7366b89, 0x3370d13e31b2c8cb}}, {{0xf768e4ca817038c3, 0x0212e5ec52afd79c, 0x4e9b401b69512c25, 0x963351c9c7b6f26a, 0x7bd7fcefd237ac19, 0xf8029ae0ea845257, 0x7f29e49f8645aa68, 0x4a79ce0d3c033d70}}};

    // real
    error |= xISOG(&A, &P, &Pd, &K, 359, 0);

    fp_cswap(&P.x, &Pd.x, 1);
    fp_cswap(&P.z, &Pd.z, 1);

    uint_c cof = {{1}};
    uint_c l;
    xMUL(&K, &A, &P, &cof);
    uint_set(&l, 353);
    xMUL(&Pd, &A, &Pd, &l);

    // real
    error |= xISOG(&A, &P, &Pd, &K, 353, 0);

    fp_inv(&A.z);
    fp_mul2(&A.x, &A.z);

#endif

    char attackCheck = 5;

    if (memcmp(&A.x, &shared_secret, sizeof(public_key)))
    {
        // not equal
        attackCheck = 0;
    }
    else
    {
        // equal
        attackCheck = 1;
    }
    simpleserial_put('r', 1, (uint8_t *)&attackCheck);

    return attackCheck;
}

int main(void)
{
    platform_init();
    init_uart();
    trigger_setup();

    simpleserial_init();

    simpleserial_addcmd('g', 0, csidh_glitch);

    while (1)
        simpleserial_get();
}
