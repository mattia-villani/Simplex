# Author: Mattia Villani
# mattiavillani94@gmail.com

# followed
# http://math.uww.edu/~mcfarlat/s-prob.htm
# and
# https://www.utdallas.edu/~scniu/OPRE-6201/documents/LP06-Simplex-Tableau.pdf
# and
# http://www.uobabylon.edu.iq/eprints/publication_9_23116_31.pdf

waitPress = False
max_iteration=30

def isclose( a, b ):
    r = 0.0001
    if b == 0 or a == 0:
        return abs(a)<r and abs(b)<r
    else :
        return abs(float(a)/float(b)-1.0) < r

def simplex( A, b, fun ):
    def dotmul( a, b ):
        if len(a) != len(b) :
            raise ValueError('Scalar multiplication can\'t be performed with vectors of different sizes')
        return [a[i]*b[i] for i in range(0,len(a))] 
    def sclmul( a, b ):
        if len(a) != len(b) :
            raise ValueError('Scalar multiplication can\'t be performed with vectors of different sizes')
        return sum( dotmul(a,b), 0 ) 

    # A * x <= b and min fun * x
    n = len(A)
    m = len(A[0])
    if not all (m==len(row) for row in A):
        print("Error, the A matrix does not have constant number of columns");
        return None;
    print("The system is "+str(n)+"x"+str(m)+".")
    if len(fun) != m :
        print("Error, the fun application has different size ("+str(len(fun))+")")
        return None
    if len(b) != n:
        print("Error, the b vector has different size ("+str(len(b))+")")
        return None
    print("Creating "+str(m)+" main variable, and "+str(n)+" auxiliar ");
    var_list = []
    for i in range(1, n+m+1):
        var_list.append( ("x" if i<=m else "s")+str(i) )
    print( "Variables: x="+str(var_list) )

    base = var_list[m: n+m]
    body = [ A[i]+[ int(i==j) for j in range(0,n) ] for i in range(0,n) ]
    fun = fun+[0 for i in range(0,n)]
    obj = [1 for i in range(0,m)]+[0 for i in range(0,n)]
    obj = dotmul( fun, obj )

    img = 0#sclmul( fun, obj )
    counter = 0

    def printStep():
        global waitPress
        print "Step "+str(counter)
        print " Base|",
        for var in var_list:
            print "{v:>5s}".format(v=var),
        print "|{v:>5s}".format(v='b')
        for i in range(0,len(base)):
            print "{v:>5s}|".format(v=base[i]),
            for el in body[i]:
                print "{v:5.2G}".format(v=el),
            print"|{b:5.2G}".format(b=b[i])
        print "".ljust(1+(3+n+m)*5, "~")
        print "  Obj|",
        for el in obj:
            print "{v:5.2G}".format(v=el),
        print "|{b:5.2G}".format(b=img),
        print ""
        if waitPress :
            try:
                input("Press a key to continue....")
            except:
                pass
    def stopConditionReached():
        return counter >= max_iteration or 0==sum([ int(x<0) for x in obj ]) # counter>=m or
    
    print "body: ["
    for row in body:
        print "      ",
        for e in row:
            print "{v:5.3f}".format(v=e),
        print ""
    print "      ]"
    print "fun : "+ str(fun);
    print "STARTING"
    printStep()
    while not stopConditionReached() :            
        pivot_col = 0
        for i in range(1, len(obj)):
            if obj[i]<obj[pivot_col] :
                pivot_col = i
        print"Selected as col pivot "+str(pivot_col)
        pivot_row = -1
        min_ratio = -1
        for i in range(0,n):
            if body[i][pivot_col] > 0:
                cur_ratio = float(b[i])/float(body[i][pivot_col])
                if (pivot_row == -1) or (cur_ratio <= min_ratio and cur_ratio>=0) :
                    pivot_row = i
                    min_ratio = cur_ratio
        if pivot_row == -1 :
            raise ValueError('Met a column of only 0s or negative value (the column '+str(pivot_col)+')')
        else :
            print"Selected as row pivot "+str(pivot_row)
        coef = float(body[pivot_row][pivot_col])
        body[pivot_row] = [ float(x)/coef for x in body[pivot_row] ]
        b[pivot_row] = b[pivot_row] / coef 
        base[pivot_row]=var_list[pivot_col]
        obj_coef = obj[pivot_col]
        print "Performing obj-("+str(obj_coef)+")*R'"+str(pivot_row)
        obj = [ obj[i]-body[pivot_row][i]*obj_coef for i in range(0,len(obj)) ]
        img = img - obj_coef * b[pivot_row]
        for i in [x for x in range(0,n) if x != pivot_row]:
            print"Performing R"+str(i),
            print"-({w:f}/{x:f})*R{y:d}".format(w=body[i][pivot_col], x=coef, y=pivot_row ),
            ratio = float(body[i][pivot_col])
            print" // ratio = "+str(ratio)

            for j in range(0,len(body[i])):
                body[i][j] = body[i][j]-ratio*body[pivot_row][j]
            b[i] = b[i] - ratio*b[pivot_row]

        obj = [ e if not isclose(e,0) else 0 for e in obj ]
        b = [ e if not isclose(e,0) else 0 for e in b ]
        body = [ [ e if not isclose(e,0) else 0 for e in row ] for row in body ]
        
        counter = counter+1
        printStep()

    return [img, [ b[base.index(v)] if v in base else 0 for v in var_list ],obj]

def gameTheory( S ):
    print "S =   ",
    for j in range(0, len(S[0])):
        print "{s:>5s}".format(s="y"+str(j+1)),
    print ""
    for i in range(0,len(S)):
        print "{s:>5s}|".format(s="x"+str(i+1)),
        for j in range(0, len(S[i])):
                print "{v:5.2f}".format(v=S[i][j]),
        print "|"
    print "BPlayer(Y): Max 1/V = ",
    for i in range(0, len(Y)):
        print "y"+str(1+i),
    print "where V is Min and y<i> = p<i>/V"
    fun = [ -1 for i in range (0,len(S[0])) ]
    print "Target fun to min "+str(fun)
    b = [1 for i in range(0,len(S))]
    print "B vector "+str(b)
    Am = S#+[ [ -int(i==j) for j in range(0,len(Y)) ] for i in range(0,len(Y)) ]
    result = simplex(Am,b,fun)
    img = result[0]
    var_res=result[1]
    obj_res=result[2]
    print "Got result of game as " 
    print "yi = pi/V = "+str(["{v:5.4f}".format(v=r) for r in var_res[:len(S[0])]])
    print "xi = qi/V = "+str(["{v:5.4f}".format(v=r) for r in obj_res[len(S[0]):]])
    print "Where 1/V = "+str(result[0])
    V = 1/img
    print "  =>  V = "+str(V)
    p =  [ V * r for r in var_res[:len(S[0])]]
    q =  [ V * r for r in obj_res[len(S[0]):]]
    print "pi = "+str(["{v:5.4f}".format(v=r) for r in p])+" // sum = {s:5.4f}".format(s=sum(p))
    print "qi = "+str(["{v:5.4f}".format(v=r) for r in q])+" // sum = {s:5.4f}".format(s=sum(q))
    maxp = max( p )
    maxq = max( q )
    print "Max p "+str(maxp)
    print "Max q "+str(maxq)
    return ((maxq,q), (maxp,p))
    

#Tests
Problems = [ # name, A , b, fun to min, expectedResult(with format [V|values])
        [
            None,
            [[2,3], [-3,2], [0,2], [2,1]],
            [ 6, 3, 5, 4 ],
            [ -4, -3 ],
            [9.0, [3.0/2.0, 1, 0, 11.0/2.0, 3.0, 0],[0, 0, 0.5, 0, 0, 1.5]]
        ],
        [
            None,
            [[2,1,1], [4,2,3], [2,5,5]],
            [ 14, 28, 30 ],
            [ -1, -2, 1 ],
            [13, [5, 4, 0, 0, 0, 0], [0, 0, 3.0625, 0, 0.0625, 0.375]]
        ],
        [
            None,
            [[-1,1], [1,-2]],
            [ 1, 2 ],
            [ 2, 1 ],
            [0, [0, 0, 1, 2], [2, 1, 0, 0] ]
        ],
        [
            "Example 1",
            [[3,-2,4], [-1,4,2], [2,2,6]],
            [1,1,1],
            [-1,-1,-1],
            [0.5, [0.4, 0.1, 0, 0, 1.0, 0], [0, 0, 2.0, 0, 0, 0.5]]
        ]
    ]

#Testing

cc = 1
failed = 0
unverified = 0
for P in Problems :        
    print ""
    print ""
    print "TEST "+str(cc)+(P[0] if P[0] is not None else "")
    cc = cc+1
    A = P[1]
    b = P[2]
    fun = P[3]
    expRes = P[4]
    result = simplex(A,b,fun)
    print "Result ( Z = "+str(result[0])+" ): "
    for i in range(0,len(result[1])):        
        print "     "+("x" if i<len(fun) else "s")+str(i)+" = "+str(result[1][i])
    print "     objs = "+str(result[2])  
    if expRes is not None:
        print "Expected result was "+str(expRes)+".... ", 
        if len(expRes) != len(result) or len(expRes[1]) != len(result[1]) or len(expRes[2]) != len(result[2]):
            print "FAILED (badly formatted result value)"
            failed = failed + 1
        elif not isclose(result[0], expRes[0] ):
            print "FAILED (unmatching Z value)"
            failed = failed + 1
        elif sum( [ int(not isclose(result[1][i],expRes[1][i])) for i in range(0, len(expRes[1]) ) ] )!= 0 :
            print "FAILED (unmatching x's values)"
            failed = failed + 1
        elif sum( [ int(not isclose(result[2][i],expRes[2][i])) for i in range(0, len(expRes[2]) ) ] )!= 0 :
            print "FAILED (unmatching obj's values)"
            failed = failed + 1
        else:
            print "OK"
    else :
        unverified = 1 + unverified

waitPress = True
print("---End. Tests failed "+str(failed)+"; unverified "+str(unverified)+"; total "+str(cc-1)+" ---")
if failed != 0 :
    print "Exiting because failure occured"
    exit()
elif waitPress:
    try:
        input("Press a key to start the colonel blotto game")
    except:
        pass


print ""
print "COLONEL BLOTTO"
print ""
xNum = 4
yNum = 3
X = [ (x,xNum-x) for x in range(0,xNum+1) ]
Y = [ (y,yNum-y) for y in range(0,yNum+1) ]
print "strategy X="+str(X)
print "strategy Y="+str(Y)
def stVal(x1,y1, x2, y2):
    def stField(x,y):
        if x>y:
            return 1+y
        elif y>x:
            return -1-x
        else:
            return 0
    return stField(x1,y1) + stField(x2, y2)
S = [ [ stVal(xx[0], yy[0], xx[1], yy[1]) for yy in Y ] for xx in X ]
((maxq,q), (maxp,p)) = gameTheory(S)
print "Blotto choses "+str(X[q.index(maxq)])+" because prob "+str(maxq)
print "Kyle   choses "+str(Y[p.index(maxp)])+" because prob "+str(maxp)


if waitPress:
    try:
        input("Press a key to start the colonel blotto reduced game")
    except:
        pass

print ""
print "COLONEL BLOTTO Reduced"
print ""
xNum = 4
yNum = 3
X = [ (x,xNum-x) for x in range(0,xNum+1) ]
Y = [ (y,yNum-y) for y in range(0,yNum+1) ]
print "strategy X="+str(X)
print "strategy Y="+str(Y)
def stVal(x1,y1, x2, y2):
    def stField(x,y):
        if x>y:
            return 1+y
        elif y>x:
            return -1-x
        else:
            return 0
    return stField(x1,y1) + stField(x2, y2)
S = [ [ stVal(xx[0], yy[0], xx[1], yy[1]) for yy in Y ] for xx in X ]
rs = len(S)
cs = len(S[0])
S1 = [ [ (S[i][j]+S[rs-i-1][j]+S[i][cs-j-1]+S[rs-i-1][cs-j-1])/4.0 for j in range(0,(cs+1)/2)] for i in range(0,(rs+1)/2)]
((maxq,q), (maxp,p)) = gameTheory(S1)
print "Blotto choses "+str(X[q.index(maxq)])+" or "+str(X[len(X)-1-q.index(maxq)])+" because prob "+str(maxq)+" // each with prob "+str(maxq/2.0)
print "Kyle   choses "+str(Y[p.index(maxp)])+" or "+str(Y[len(Y)-1-p.index(maxp)])+" because prob "+str(maxp)+" // each with prob "+str(maxq/2.0)


