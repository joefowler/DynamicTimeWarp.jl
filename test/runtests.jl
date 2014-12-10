using DynamicTimeWarp
using Base.Test

# write your own tests here
@test 1 == 1


# Test the Distance functions
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


# Test dtw itself
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

