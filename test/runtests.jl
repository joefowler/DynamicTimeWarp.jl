using DynamicTimeWarp
using Base.Test

#############################################
# Test the Distance functions
#############################################

Distance = DynamicTimeWarp.Distance
@test Distance.square(0,0) == 0
@test Distance.square(0,1) == 1
@test Distance.square(9,10) == 1
@test Distance.square(0,5) == 25
@test Distance.square(0.0,5.0) == 25.0
@test Distance.square(0.0,5) == 25.0

@test Distance.absval(0,0) == 0
@test Distance.absval(0,3) == 3
@test Distance.absval(0,5) == 5
@test Distance.absval(-2,1) == 3
@test Distance.absval(-2,3) == 5
@test Distance.absval(-2,5) == 7
@test Distance.absval(-2,-5) == 3

s1, s2 = [1:5], [2:2:10]
pc = Distance.poissonclosure(s1, s2)
@test pc(1,2) == 0
@test pc(5,10) == 0
@test pc(5,5) == (5/sum(s1)-5/sum(s2))^2
@test pc(5,0) == (5/sum(s1))^2
@test pc(0,5) == (5/sum(s2))^2


#############################################
# Test dtw itself
#############################################

a=[1,1,1,2,4,6,5,5,5,4,4,3,1,1,1]
b=[1,1,2,4,6,6,6,5,4,4,4,3,3,3,1]
cost, match1, match2 = dtw(a,b)
@test cost==0
@test match1==[1,2,3,4,5,6,6,6,7,8,9,10,10,11,12,12,12,13,14,15]
@test match2==[1,1,2,3,4,5,6,7,8,8,8, 9,10,11,12,13,14,15,15,15]

a[end] += 2
cost, match1, match2 = dtw(a,b)
@test cost==4


a=[1:10]
b=a+1
cost, match1, match2 = dtw(a,b)
@test cost==2

a=zeros(Int,6)
b=1+a
cost, match1, match2 = dtw(a,b)
@test cost==length(a)

# Verify that a tie prefers diagonal moves
a=[1,1,1]
b=[1,1,1]
cost, pa, pb = dtw(a,b)
@test cost==0
@test pa==[1,2,3]
@test pb==[1,2,3]

# Verify that trackback ends properly if it reaches an edge before reaching [1,1]
# Also check that trackback prefers diagonal moves
a=[0,1,1,1]
b=[0,0,1,1]
cost, pa, pb = dtw(a,b)
@test cost==0
@test pa==[1,1,2,3,4]
@test pb==[1,2,3,3,4]



#############################################
# Test DTW with windows
#############################################

# Verify that a tie prefers diagonal moves
a=[1,1,1]
b=[1,1,1]
cost, pa, pb = dtwwindowed(a,b,[1,1,1],[3,3,3])
@test cost==0
@test pa==[1,2,3]
@test pb==[1,2,3]

# Verify that trackback ends properly if it reaches an edge before reaching [1,1]
# Also check that trackback prefers diagonal moves
a=[0,1,1,1]
b=[0,0,1,1]
cost, pa, pb = dtwwindowed(a,b,[1,1,1,1],[4,4,4,4])
@test cost==0
@test pa==[1,1,2,3,4]
@test pb==[1,2,3,3,4]

# First do the windowed test w/o windows
a=[0,1,2,3,4,4,4,4]
b=[0,0,1,2,2,2,3,4]
best_pa = [1,1,2,3,3,3,4,5,6,7,8]
best_pb = [1,2,3,4,5,6,7,8,8,8,8]
cost, pa, pb = dtw(a,b)
@test cost == 0
@test pa == best_pa
@test pb == best_pb

# Wide window, not touching optimal path
rmin = [1,1,1,2,3,4,5,6]
rmax = [4,6,7,8,8,8,8,8]
cost, pa, pb = dtwwindowed(a,b,rmin,rmax)
@test cost == 0
@test pa == best_pa
@test pb == best_pb

# Bottom of window is optimal path
rmin = [1,3,4,7,8,8,8,8]
rmax = [4,6,7,8,8,8,8,8]
cost, pa, pb = dtwwindowed(a,b,rmin,rmax)
@test cost == 0
@test pa == best_pa
@test pb == best_pb

# Top of window is optimal path
rmin = [1,1,1,2,3,4,5,6]
rmax = [2,3,6,7,8,8,8,8]
cost, pa, pb = dtwwindowed(a,b,rmin,rmax)
@test cost == 0
@test pa == best_pa
@test pb == best_pb

# Top and bottom of window are optimal path
rmin = [1,3,4,7,8,8,8,8]
rmax = [2,3,6,7,8,8,8,8]
cost, pa, pb = dtwwindowed(a,b,rmin,rmax)
@test cost == 0
@test pa == best_pa
@test pb == best_pb

# Now top of window cuts into optimal path
rmin = [1,1,1,2,3,4,5,6]
rmax = [4,4,5,6,7,8,8,8]
cost, pa, pb = dtwwindowed(a,b,rmin,rmax)
@test cost == 2
@test pa == [1,1,2,3,3,4,5,6,7,8]
@test pb == [1,2,3,4,5,6,7,8,8,8]


#############################################
# Test sequence compressions
#############################################

compress = DynamicTimeWarp.compress
s=[0:2:98]
s1 = compress(s)
s2 = compress(s1)
@test s1==float([1:4:97])
@test s2==float(vcat([3:8:91],[97]))

s = [1]
s1 = compress(s)
@test s1==[1.0]


#############################################
# Test window computation
#############################################
computewindow = DynamicTimeWarp.computewindow

# Simplest path (along the diagonal)
p=[1:8]
rmin,rmax = computewindow(p,p,1)
@test rmin==[1,1,1,2,3,4,5,6]
@test rmax==[3,4,5,6,7,8,8,8]

rmin,rmax = computewindow(p,p,2)
@test rmin==[1,1,1,1,1,2,3,4]
@test rmax==[5,6,7,8,8,8,8,8]

# A warpy path
pa=[1,1,2,3,4,5,6,7,8,8,8]
pb=[1,2,3,3,3,4,4,5,6,7,8]
rmin,rmax = computewindow(pa, pb, 1)
@test pa[end]==length(rmin)
@test pa[end]==length(rmax)
@test rmin==[1,1,2,2,2,3,3,4]
@test rmax==[4,4,4,5,5,6,8,8]

rmin,rmax = computewindow(pa, pb, 2)
@test pa[end]==length(rmin)
@test pa[end]==length(rmax)
@test rmin==[1,1,1,1,1,1,2,2]
@test rmax==[5,5,6,6,7,8,8,8]

rmin,rmax = computewindow(pa, pb, 20)
@test pa[end]==length(rmin)
@test pa[end]==length(rmax)
@test rmin==fill(1,8)
@test rmax==fill(8,8)

# Extreme path: follows left then upper edge
pa=[1,1,1,1,1,1,1,1,2,3,4,5,6,7,8]
pb=[1,2,3,4,5,6,7,8,8,8,8,8,8,8,8]
rmin,rmax = computewindow(pa, pb, 1)
@test pa[end]==length(rmin)
@test pa[end]==length(rmax)
@test rmin==[1,1,7,7,7,7,7,7]
@test rmax==[8,8,8,8,8,8,8,8]

rmin,rmax = computewindow(pa, pb, 2)
@test pa[end]==length(rmin)
@test pa[end]==length(rmax)
@test rmin==[1,1,1,6,6,6,6,6]
@test rmax==[8,8,8,8,8,8,8,8]


# More columns than rows
pa=[1,2,3,4,5,6,7,8]
pb=[1,2,3,4,4,4,4,4]
rmin,rmax = computewindow(pa, pb, 1)
@test rmin==[1,1,1,2,3,3,3,3]
@test rmax==[3,4,4,4,4,4,4,4]

rmin,rmax = computewindow(pa, pb, 2)
@test rmin==[1,1,1,1,1,2,2,2]
@test rmax==fill(4,8)

rmin,rmax = computewindow(pa, pb, 3)
@test rmin==fill(1,8)
@test rmax==fill(4,8)

rmin,rmax = computewindow(pa, pb, 4)
@test rmin==fill(1,8)
@test rmax==fill(4,8)

rmin,rmax = computewindow(pa, pb, 47)
@test rmin==fill(1,8)
@test rmax==fill(4,8)


# More rows than columns
pa,pb = pb,pa
rmin,rmax = computewindow(pa, pb, 1)
@test rmin==[1,1,1,2]
@test rmax==[3,4,8,8]

rmin,rmax = computewindow(pa, pb, 2)
@test rmin==fill(1,4)
@test rmax==[5,8,8,8]

rmin,rmax = computewindow(pa, pb, 3)
@test rmin==fill(1,4)
@test rmax==fill(8,4)

rmin,rmax = computewindow(pa, pb, 4)
@test rmin==fill(1,4)
@test rmax==fill(8,4)

rmin,rmax = computewindow(pa, pb, 47)
@test rmin==fill(1,4)
@test rmax==fill(8,4)

