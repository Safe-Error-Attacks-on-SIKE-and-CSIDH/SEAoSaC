# #     bit = (m[9 >> LOG2RADIX] >> (9 & (RADIX-1))) & 1;
    
# #ifdef ATTACK_BIT_9
#     trigger_high(); 
#     swap = bit ^ prevbit;
#     prevbit = bit;
#     mask = 0 - (digit_t)swap;
#     trigger_low(); 
# #else
#     swap = bit ^ prevbit;
#     prevbit = bit;
#     mask = 0 - (digit_t)swap;
# #endif 

# #     swap_points(R, R2, mask);
# #     xDBLADD(R0, R2, R->X, A24);        
# #     fp2mul_mont(R2->X, R->Z, R2->X); 

S = ""
for i in range(0, 9):
    S = S + "\n//------------------------------------\nbit = (m[" + str(i) + " >> LOG2RADIX] >> (" + str(i) + " & (RADIX-1))) & 1;\n\n#ifdef ATTACK_BIT_" + str(i) + "\n    trigger_high();\n    swap = bit ^ revbit;\n    prevbit = bit;\n    mask = 0 - (digit_t)swap;\n    trigger_low();\n#else\n    swap = bit ^ prevbit;\n    prevbit = bit;\n    mask = 0 - (digit_t)swap;\n#endif\n\nswap_points(R, R2, mask);\nxDBLADD(R0, R2, R->X, A24);\nfp2mul_mont(R2->X, R->Z, R2->X);\n"

print(S)