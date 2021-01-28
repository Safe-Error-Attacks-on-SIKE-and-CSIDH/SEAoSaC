SCOPETYPE = 'OPENADC'
PLATFORM = 'CWLITEARM'
sample_size = 5

def sort_glitch(glitch):
    return glitch[2]

import chipwhisperer as cw
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



import time
time.sleep(0.05)
scope.default_setup()
def reset_target(scope):
    if PLATFORM == "CW303" or PLATFORM == "CWLITEXMEGA":
        scope.io.pdic = 'low'
        time.sleep(0.05)
        scope.io.pdic = 'high_z' #XMEGA doesn't like pdic driven high
        time.sleep(0.05)
    else:  
        scope.io.nrst = 'low'
        time.sleep(0.05)
        scope.io.nrst = 'high'
        time.sleep(0.05)

fw_path = "/home/sopmac/opt/chipwhisperer-5.1.3/hardware/victims/firmware/glitch-simple/glitchsimple-{}.hex".format(PLATFORM)
#fw_path = "/home/sopmac/Desktop/phd/git@github/m_safe_error/attack/speed-{}.hex".format(PLATFORM)


# cd ../hardware/victims/firmware/glitch-simple
# make PLATFORM=$1 CRYPTO_TARGET=NONE FUNC_SEL=GLITCH_INF

cw.program_target(scope, prog, fw_path)

reset_target(scope)
target.flush()
time.sleep(1)
resp = target.read()
print(resp)

from collections import namedtuple
Range = namedtuple('Range', ['min', 'max', 'step'])
if PLATFORM == "CWLITEARM" or PLATFORM == "CW308_STM32F3":
    scope.glitch.clk_src = "clkgen"
    scope.glitch.output = "glitch_only"
    scope.glitch.trigger_src = "ext_single"
    scope.glitch.width = 35
    scope.glitch.offset = -17.4
    scope.glitch.repeat = 1

    width_range = Range(38.5, 39.1, 0.4)
    offset_range = Range(-28.4, -28.125, 0.4)
    scope.glitch.offset_fine = 24
    def glitch_on(scope):
        scope.io.glitch_lp = False
        scope.io.glitch_hp = True
    def glitch_off(scope):
        scope.gio.glitch_hp = False
    glitch_on(scope)
    scope.glitch.ext_offset = 2186
    print(scope.glitch)
elif PLATFORM == "CWNANO" and SCOPETYPE == "CWNANO":
    scope.glitch.ext_offset = 546
    scope.adc.clk_freq = 7.5E6
    scope.glitch.repeat = 6
    repeat_range = range(4, 7)
    offset_range = range(475, 510)
    def glitch_on(scope):
        pass
    def glitch_off(scope):
        pass
    pass #later

print(scope)

from tqdm.notebook import tnrange
reset_target(scope)
target.flush()

# if SCOPETYPE == "OPENADC":
#     scope.glitch.trigger_src = "ext_continuous"

# for j in range(20):
#     line = ""

#     while "\n" not in line:
#         time.sleep(0.1)
#         line += target.read()
#     print(line)
#     lines = line.split("\n")
#     if len(lines) > 1:
#         line = lines[-1]
#     else:
#         line = ""

#     while "\n" not in line:
#         time.sleep(0.1)
#         line += target.read()
#     print(line)
#     if "hello" in line:
#         print("Target crashed")
#     nums = line.split(" ")
#     #print(line)
#     try:
#         if int(nums[0]) != 40000:
#             print("+" + line)
#         print("-" + line)
#     except ValueError as e:
#         continue

# if SCOPETYPE == "OPENADC":
#     scope.glitch.trigger_src = "ext_single"

from tqdm import tnrange, tqdm_notebook
reset_target(scope)
glitches = []
glitch_text = []

if SCOPETYPE == "OPENADC":
    target.flush()
    scope.glitch.trigger_src = "ext_continuous"
    scope.glitch.offset_fine = 24
    scope.glitch.repeat = 1
    scope.glitch.ext_offset = 2186
    scope.glitch.offset = offset_range.min
    t_offset = tqdm_notebook(total=int((offset_range.max-offset_range.min)/offset_range.step) + 1, desc="Offset")
    while scope.glitch.offset < offset_range.max:
        scope.glitch.width = width_range.min
        t_width = tqdm_notebook(total=int((width_range.max-width_range.min)/width_range.step), leave=False, desc="Width")
        while scope.glitch.width < width_range.max:
            successes = 0
            crashes = 0
            for j in tnrange(sample_size, leave=False, desc="Attempt"):
                line = ""
                while "\n" not in line:
                    time.sleep(0.1)
                    num_char = target.in_waiting()
                    if num_char == 0:
                        glitch_off(scope)
                        time.sleep(0.01)
                        glitch_on(scope)
                        break
                    line += target.read()
                lines = line.split("\n")
                if len(lines) > 1:
                    line = lines[-1]
                else:
                    line = ""

                while "\n" not in line:
                    time.sleep(0.1)
                    num_char = target.in_waiting()
                    if num_char == 0:
                        glitch_off(scope)
                        time.sleep(0.01)
                        glitch_on(scope)
                        break
                    line += target.read()

                nums = line.split(" ")
                if "hello" in line:
                    crashes += 1
                    #print("Target crashed")
                #print(line)
                try:
                    if nums[0] == "":
                        continue
                    if int(nums[0]) != 40000:
                        glitch_text += line
                        successes += 1
                except ValueError as e:
                    continue
            glitches.append([scope.glitch.width, scope.glitch.offset, successes / sample_size, crashes / sample_size])
            if successes > 0:
                print([scope.glitch.width, scope.glitch.offset, successes / sample_size, crashes / sample_size])
            scope.glitch.width += width_range.step
            t_width.update()

        scope.glitch.offset += offset_range.step
        t_offset.update()
        t_width.close()
    t_width.close()
    t_offset.close()
    scope.glitch.trigger_src = "ext_single"    

glitches.sort(key=sort_glitch,reverse=True)
for glitch in glitches:
    print(glitch)
