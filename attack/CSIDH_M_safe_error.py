import os
import datetime
import sys
import chipwhisperer as cw
import time
import chipwhisperer.common.results.glitch as glitch
from importlib import reload
import chipwhisperer.common.results.glitch as glitch
from tqdm.notebook import trange
import struct

def main(argv):
    if(len(argv) != 4):
        print(
            "\nplease specify the path, attack dummy=0 or real=1 isogeny, random=0 or critical spot=1, and numer of trials\npython CSIDH_M_safe_error_random.py [path] [dummy/real] [random/spot] [ntrials]\ne.g. python M_safe_error.py /home/me/chipwhisperer/hardware/victims/firmware/csidh/ 0 1 10\n")
        sys.exit()
    else:
        PATH = argv[0]
        KEY = argv[1]
        TYPE = int(argv[2])
        NTRIAL = int(argv[3])

        make = "make PLATFORM=CWLITEARM CRYPTO_TARGET=NONE FUNC_SEL=KEY_" + KEY        

        os.chdir(PATH)
        os.system(make)

        SCOPETYPE = 'OPENADC'
        PLATFORM = 'CWLITEARM'

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

        fw_path = PATH + "csidh512-{}.hex".format(PLATFORM)

        cw.program_target(scope, prog, fw_path)

        scope.clock.clkgen_freq = 24E6
        target.baud = 38400*24/7.37
        time.sleep(0.1)
        def reboot_flush():
            scope.io.nrst = False
            time.sleep(0.05)
            scope.io.nrst = "high_z"
            time.sleep(0.05)
            #Flush garbage too
            target.flush()

        reboot_flush()
        scope.arm()
        target.write("g\n")

        scope.capture()
        time.sleep(7)
        val = target.simpleserial_read_witherrors('r', 1, glitch_timeout=10)#For loop check
        valid = val['valid']
        if valid:
            response = val['payload']
            raw_serial = val['full_response']
            error_code = val['rv']
            if val['rv'] == 1:
                if TYPE == 0:
                    print("All checks passed.😍\nAttacking randomly ...")
                else:
                    print("All checks passed.😍\nAttacking critical spot ...")
            else:
                print("Communication and/or key check error!\n")
                sys.exit()
        #print(val)

        gc = glitch.GlitchController(groups=["success", "reset", "normal"], parameters=["width", "offset", "ext_offset"])

        scope.glitch.clk_src = "clkgen" # set glitch input clock
        scope.glitch.output = "glitch_only" # glitch_out = clk ^ glitch
        scope.glitch.trigger_src = "ext_single" # glitch only after scope.arm() called
        if PLATFORM == "CWLITEXMEGA":
            scope.io.glitch_lp = True
            scope.io.glitch_hp = True
        elif PLATFORM == "CWLITEARM":
            scope.io.glitch_lp = True
            scope.io.glitch_hp = True
        elif PLATFORM == "CW308_STM32F3":
            scope.io.glitch_hp = True
            scope.io.glitch_lp = True

        g_step = 0.4
        if PLATFORM=="CWLITEXMEGA":
            gc.set_range("width", 45.7, 47.8)
            gc.set_range("offset", 2.8, 10)
            gc.set_range("ext_offset", 2, 4)
            scope.glitch.repeat = 10
        elif PLATFORM == "CWLITEARM":
            #should also work for the bootloader memory dump
            gc.set_range("width", 34.7, 36)
            gc.set_range("offset", -41, -30)
            gc.set_range("ext_offset", 0, 40)
            scope.glitch.repeat = 7
        elif PLATFORM == "CW308_STM32F3":
            #these specific settings seem to work well for some reason
            #also works for the bootloader memory dump
            gc.set_range("ext_offset", 9, 12)
            gc.set_range("width", 47.6, 49.6)
            gc.set_range("offset", -19, -21.5)
            scope.glitch.repeat = 5

        gc.set_global_step(g_step)

        scope.adc.timeout = 5

        glitch_cnt = 0

        reboot_flush()
        sample_size = 1
        loff = scope.glitch.offset
        lwid = scope.glitch.width
        total_successes = 0

        if TYPE == 0:    # RANDOM
            for glitch_setting in gc.glitch_values():            
                print("{},".format(glitch_cnt), end='')
                sys.stdout.flush()
                scope.glitch.offset = glitch_setting[1]
                scope.glitch.width = glitch_setting[0]
                scope.glitch.ext_offset = glitch_setting[2]
         
                successes = 0
                resets = 0
                for i in range(NTRIAL):
                    glitch_cnt = glitch_cnt + 1
                    reboot_flush()
                    target.flush()
                    if scope.adc.state:
                        reboot_flush()
                        resets += 1

                    scope.arm()

                    target.write("g\n")

                    ret = scope.capture()
                    time.sleep(7)

                    val = target.simpleserial_read_witherrors('r', 1, glitch_timeout=10)#For loop check

                    scope.io.glitch_hp = False
                    scope.io.glitch_hp = True
                    scope.io.glitch_lp = False
                    scope.io.glitch_lp = True
                    if ret:
                        resets += 1
                        reboot_flush()

                    else:
                        if val['valid'] is False:
                            reboot_flush()
                            resets += 1
                        else:
                            if val['rv'] == 0: #for loop check
                                successes += 1
                if successes > 0:
                    print("successes = {}, resets = {}, offset = {}, width = {}, ext_offset = {}".format(successes, resets, scope.glitch.offset, scope.glitch.width, scope.glitch.ext_offset))
                    total_successes += successes

        else:
            for _ in range(100):
                print("{},".format(glitch_cnt), end='')
                sys.stdout.flush()

                # critical spot
                scope.glitch.offset = -33.59375
                scope.glitch.width = 35.546875
                scope.glitch.ext_offset = 2

                successes = 0
                resets = 0
                for i in range(NTRIAL):
                    glitch_cnt = glitch_cnt + 1
                    reboot_flush()
                    target.flush()
                    if scope.adc.state:
                        reboot_flush()
                        resets += 1

                    scope.arm()
                    target.write("g\n")

                    ret = scope.capture()
                    time.sleep(7)

                    val = target.simpleserial_read_witherrors('r', 1, glitch_timeout=10)#For loop check

                    #print(val)
                    scope.io.glitch_hp = False
                    scope.io.glitch_hp = True
                    scope.io.glitch_lp = False
                    scope.io.glitch_lp = True
                    if ret:
                        resets += 1
                        reboot_flush()

                    else:
                        if val['valid'] is False:
                            reboot_flush()
                            resets += 1
                        else:
                            if val['rv'] == 0: #for loop check
                                successes += 1
                if successes > 0:
                    total_successes += successes

        print("total successes = {}, trials = {}, ratio = {}".format(total_successes, glitch_cnt, total_successes/glitch_cnt))

        scope.dis()
        target.dis()

if __name__ == "__main__":
    main(sys.argv[1:])