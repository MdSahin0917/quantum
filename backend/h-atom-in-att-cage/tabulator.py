# v0_values = [-1.0,-0.9,-0.8,-0.7,-0.6,-0.55,-0.5,-0.48,-0.45,-0.43,-0.4,-0.3,-0.2,-0.1,-0.08,-0.05]
#v0_values = [-1.0,-0.8,-0.6,-0.4,-0.2,-0.125,-0.08,-0.05,-0.01]
#v0_values = [-0.001,-0.005,-0.008,-0.01,-0.05,-0.055556,-0.08,-0.1]
import sys
dict = {
  0: "s",
  1: "p",
  2: "d",
  3: "f",
  4: "g",
  5: "h"
}

if len(sys.argv) < 5:
    print("Usage: python child.py v0 n l N")
    sys.exit(1)

v0 = float(sys.argv[1])
n = int(sys.argv[2])
l = int(sys.argv[3])
N = int(sys.argv[4])
maxR = int(sys.argv[5])

print(f"Running tabulator with v0 = {v0}, n = {n}, l = {l}, N = {N}")




# for v0 in v0_values:
#     File = f"{n}{dict[l]}_v0{v0}.dat"

   


#     loop = [0.1 + 0.1 *i for i in range(0, 301)] #[0.5 + 0.1 * i for i in range(0, 101)]
#     fi = open(File,'w')
#     # fi.write('x0 ,  EN   , KIN , POT,   COUL, CENT ,  r , r^2 , r^3  , 1/r , 1/r^2 , 1/r^3 ,SD(r) ,PC(r), kirk, buck, sq, tun,swell,Sr, Dr, LMC\n')
#     for x0 in loop:
#         V0 = v0
#         D0 = 5.0
#         x0 = round(x0,3)
#         name_folder = f'(V = {V0}, R = {x0} , Delta = {D0})'
#         f = open(f'{n}{dict[l]}_v0{v0}/{name_folder}/Correlated_data_(V = {V0}, R = {x0} , Delta = {D0}).txt','r')
#         text_store = []
#         for i in f.readlines():
#             text_store.append(i)
#         f.close()
#         BigString = f'{x0}  '
#         for j in range(20,70):
#             string = text_store[j-1].split()[-1]
#             BigString += string + '  '
#         fi.write(BigString + "\n")
#     fi.close()



File = f"{n}{dict[l]}_v0{v0}.dat"

   


loop = [0.1 + 0.1 *i for i in range(0, maxR)] #[0.5 + 0.1 * i for i in range(0, 101)]
fi = open(File,'w')
# fi.write('x0 ,  EN   , KIN , POT,   COUL, CENT ,  r , r^2 , r^3  , 1/r , 1/r^2 , 1/r^3 ,SD(r) ,PC(r), kirk, buck, sq, tun,swell,Sr, Dr, LMC\n')
for x0 in loop:
    V0 = v0
    D0 = 5.0
    x0 = round(x0,3)
    name_folder = f'(V = {V0}, R = {x0} , Delta = {D0})'
    f = open(f'{n}{dict[l]}/{name_folder}/Correlated_data_(V = {V0}, R = {x0} , Delta = {D0}).txt','r')
    text_store = []
    for i in f.readlines():
        text_store.append(i)
    f.close()
    BigString = f'{x0}  '
    for j in range(20,70):
        string = text_store[j-1].split()[-1]
        BigString += string + '  '
    fi.write(BigString + "\n")
fi.close()

