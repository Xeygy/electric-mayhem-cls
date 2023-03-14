import json
import time
from scipy import stats
from multiprocessing import Process, Queue

f = open('stm32f0_aes.json')
traces = json.load(f)

plaintexts = []
ciphertexts = []
powertraces = []
key = []

# Rijndael S-box (AES's S-box)
sbox = [
    0x63, 0x7C, 0x77, 0x7B, 0xF2, 0x6B, 0x6F, 0xC5, 0x30, 0x01, 0x67, 0x2B, 0xFE, 0xD7,
    0xAB, 0x76, 0xCA, 0x82, 0xC9, 0x7D, 0xFA, 0x59, 0x47, 0xF0, 0xAD, 0xD4, 0xA2, 0xAF,
    0x9C, 0xA4, 0x72, 0xC0, 0xB7, 0xFD, 0x93, 0x26, 0x36, 0x3F, 0xF7, 0xCC, 0x34, 0xA5,
    0xE5, 0xF1, 0x71, 0xD8, 0x31, 0x15, 0x04, 0xC7, 0x23, 0xC3, 0x18, 0x96, 0x05, 0x9A,
    0x07, 0x12, 0x80, 0xE2, 0xEB, 0x27, 0xB2, 0x75, 0x09, 0x83, 0x2C, 0x1A, 0x1B, 0x6E,
    0x5A, 0xA0, 0x52, 0x3B, 0xD6, 0xB3, 0x29, 0xE3, 0x2F, 0x84, 0x53, 0xD1, 0x00, 0xED,
    0x20, 0xFC, 0xB1, 0x5B, 0x6A, 0xCB, 0xBE, 0x39, 0x4A, 0x4C, 0x58, 0xCF, 0xD0, 0xEF,
    0xAA, 0xFB, 0x43, 0x4D, 0x33, 0x85, 0x45, 0xF9, 0x02, 0x7F, 0x50, 0x3C, 0x9F, 0xA8,
    0x51, 0xA3, 0x40, 0x8F, 0x92, 0x9D, 0x38, 0xF5, 0xBC, 0xB6, 0xDA, 0x21, 0x10, 0xFF,
    0xF3, 0xD2, 0xCD, 0x0C, 0x13, 0xEC, 0x5F, 0x97, 0x44, 0x17, 0xC4, 0xA7, 0x7E, 0x3D,
    0x64, 0x5D, 0x19, 0x73, 0x60, 0x81, 0x4F, 0xDC, 0x22, 0x2A, 0x90, 0x88, 0x46, 0xEE,
    0xB8, 0x14, 0xDE, 0x5E, 0x0B, 0xDB, 0xE0, 0x32, 0x3A, 0x0A, 0x49, 0x06, 0x24, 0x5C,
    0xC2, 0xD3, 0xAC, 0x62, 0x91, 0x95, 0xE4, 0x79, 0xE7, 0xC8, 0x37, 0x6D, 0x8D, 0xD5,
    0x4E, 0xA9, 0x6C, 0x56, 0xF4, 0xEA, 0x65, 0x7A, 0xAE, 0x08, 0xBA, 0x78, 0x25, 0x2E,
    0x1C, 0xA6, 0xB4, 0xC6, 0xE8, 0xDD, 0x74, 0x1F, 0x4B, 0xBD, 0x8B, 0x8A, 0x70, 0x3E,
    0xB5, 0x66, 0x48, 0x03, 0xF6, 0x0E, 0x61, 0x35, 0x57, 0xB9, 0x86, 0xC1, 0x1D, 0x9E,
    0xE1, 0xF8, 0x98, 0x11, 0x69, 0xD9, 0x8E, 0x94, 0x9B, 0x1E, 0x87, 0xE9, 0xCE, 0x55,
    0x28, 0xDF, 0x8C, 0xA1, 0x89, 0x0D, 0xBF, 0xE6, 0x42, 0x68, 0x41, 0x99, 0x2D, 0x0F,
    0xB0, 0x54, 0xBB, 0x16,
]

# Populate lists
for trace in traces:
    plaintexts.append(trace['pt'])
    ciphertexts.append(trace['ct'])
    powertraces.append(trace['pm'])

num_traces = len(traces)
print("number of powertraces:", num_traces)  # 50

num_bytes = len(plaintexts[0])
print("Bytes in plaintext:", num_bytes) # 16

num_points = len(powertraces[0])
print("Datapoints in powertraces:", num_points) # 1806

#########################
#        HELPER         #
#       FUNCTIONS       #
#########################

# returns the number of 1's in the binary form of the input number
def hamming(input):
    return bin(input).count('1')

# returns a list of 50 values, each being the i'th point in the set 
# of 50 powertraces we've loaded in
def get_points_at_loc(i):
    res = []
    for trace in powertraces:
        res.append(trace[i])
    return res
    
# finds the best byte value for a given byte location (0-15) in the key (by correlation)
def find_best_value(loc):
    best_value = 0
    bestest_correlation = 0
    
    # try all possible byte values
    for i in range(256):
        # print out progress (approx percent)
        # if i % 25 == 0:
        #     print(f"{i//25 * 10}% done.")

        # list of 50 hamming distances for each powertrace
        hamming_array = []
        for pt in plaintexts:
            hamming_array.append(hamming(sbox[pt[loc] ^ i]))

        # best correlation stat for value i
        best_correlation = 0
        for j in range(num_points): # number of datapoints in each power trace
            points = get_points_at_loc(j)
            corr = abs(stats.pearsonr(hamming_array, points)[0])
            if best_correlation < corr:
                best_correlation = corr
        
        # update best_value if current i has the best correlation so far
        if bestest_correlation < best_correlation:
            bestest_correlation = best_correlation
            best_value = i

    return (best_value, bestest_correlation)

#########################
#          END          #
#########################


########################################################
#   Solution if you don't want to do multiprocessing   #
########################################################
'''
for i in range(num_bytes):
    print(find_best_value(i))
'''

################################
#   Multiprocessing solution   #
################################

N_PARTS = 8  # 8 processes 

# function to find the best value with a multiprocessing solution
def multiprocessing(q, start_time):
    while not q.empty():
        try:
            loc = q.get(timeout=5)
            best_value = find_best_value(loc)
            print(f"byte #{loc}, with value {best_value[0]}, correlation {best_value[1]}")
            t = (time.time() - start_time)
            print(f"time: {t//60}m {t%60}s")
        except TimeoutError:
            print("We lacked patience and got a multiprocessing.TimeoutError")

if __name__ == '__main__':
    loc_queue = Queue() # queue of 16 jobs, one for each byte of the key
    for i in range(num_bytes):
        loc_queue.put_nowait(i)
        
    procs = []
    start_time = time.time()
    for i in range(N_PARTS):
        proc = Process(target=multiprocessing, args=(loc_queue, start_time))
        procs.append(proc)
        proc.start()
    for proc in procs:
        proc.join()
