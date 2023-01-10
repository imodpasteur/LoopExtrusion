import random
print ("# initial coordinates for a polymer")
print("\n")
L=600
print("\n")
print(L, "atoms")
print (L-1, "bonds")
print (L-2, "angles")

print("\n")

print ("1  atom types")
print ("1  bond types")
print ("1  angle types")

print("\n")

print ("-22.0000   22.0000 xlo xhi")
print ("-22.0000   22.0000 ylo yhi")
print ("-22.0000   22.0000 zlo zhi")


print("\n")
print("Masses" "\n")
print("1	1.00")


print("\n")
print ("Atoms" "\n")
x=[None]*(L+1)
y=[None]*(L+1)
z=[None]*(L+1)
for i in range (1,L+1):
    x[i]=round(i/52,2)+round(random.uniform(-0.5,0.5), 2)-5.769
    y[i]=round(i/52,2)+round(random.uniform(-0.5,0.5), 2)-5.769
    z[i]=round(i/52,2)+round(random.uniform(-0.5,0.5), 2)-5.769
    print (i, 1, 1, x[i], y[i], z[i])


print("\n")
print ("Bonds" "\n")
for i in range (1,L):
    print(i, 1, i, i+1)


print("\n")
print ("Angles" "\n")
for i in range (1,L-1):
    print (i, 1, i, i+1, i+2)
