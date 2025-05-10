def cf(x,lim):
    qp=floor(x)
    exp=[qp]
    for i in range(lim):
        x=1/(x-qp)
        qp=floor(x)
        exp.append(qp)
    return exp

def convergents_cf(exp,lim):
	A=[exp[0],exp[1]*exp[0]+1]
	B=[1,exp[1]]

	for i in range(2,lim):
		A.append(exp[i]*A[i-1]+A[i-2])
		B.append(exp[i]*B[i-1]+B[i-2])
	
	return [A,B]


#compute the continued fraction of (ax+b)/(cx+d), x is a value in the first function Gosper, it is a sequence in the second function seq_Gosper


def Gosper(a,b,c,d,x,lim):
    exp=[]

    for i in range(lim):
        if c!=0 and (c+d)!=0:
            print(c)
            print(d)
            if floor(a/c)==floor((a+b)/(c+d)) and sign(c)==sign(c+d):
                qtrans=floor(a/c)
                exp.append(qtrans)
                
                [a,b,c,d]=[c,d,a-c*qtrans,b-d*qtrans]

            else:
                q=floor(x)

                [a,b,c,d]=[a*q+b,a,c*q+d,c]
                
                x=1/(x-q)
        else:
            q=floor(x)

            [a,b,c,d]=[a*q+b,a,c*q+d,c]
            
            x=1/(x-q)
    
    return exp

def seq_Gosper(a,b,c,d,x,lim):
    exp=[]

    j=0

    for i in range(lim):
        if c!=0 and (c+d)!=0:
            if floor(a/c)==floor((a+b)/(c+d)) and sign(c)==sign(c+d):
                qtrans=floor(a/c)
                exp.append(qtrans)
                
                [a,b,c,d]=[c,d,a-c*qtrans,b-d*qtrans]

            else:
                q=x[j]

                [a,b,c,d]=[a*q+b,a,c*q+d,c]
                
                j+=1
        else:
            q=x[j]

            [a,b,c,d]=[a*q+b,a,c*q+d,c]
            
            j+=1
    
    return exp

#compute the continued fraction of (axy+bx+cy+d)/(exy+fx+gy+h)


def Gosper2(a,b,c,d,e,f,g,h,x,y,lim):
    exp=[]

    [ix,iy]=[0,0]

    for i in range(lim):
        if e!=0 and (e+f)!=0 and (e+g)!=0 and (e+f+g+h)!=0:
            if floor(a/e)==floor((a+b)/(e+f))==floor((a+c)/(e+g))==floor((a+b+c+d)/(e+f+g+h)) and  sign(e)==sign(e+f)==sign(e+g):
                qtrans=floor(a/e)
                exp.append(qtrans)

                [a,b,c,d,e,f,g,h]=[e,f,g,h,a-e*qtrans,b-f*qtrans,c-g*qtrans,d-h*qtrans]
                
            else:
                if ix==iy:
                    q=floor(x)

                    [a,b,c,d,e,f,g,h]=[a*q+c,b*q+d,a,b,e*q+g,f*q+h,e,f]

                    x=1/(x-q)
                    ix+=1

                else:
                    q=floor(y)

                    [a,b,c,d,e,f,g,h]=[a*q+b,a,c*q+d,c,e*q+f,e,g*q+h,g]
            
                    y=1/(y-q)
                    iy+=1

        else:
                if ix==iy:
                    q=floor(x)

                    [a,b,c,d,e,f,g,h]=[a*q+c,b*q+d,a,b,e*q+g,f*q+h,e,f]

                    x=1/(x-q)
                    ix+=1

                else:
                    q=floor(y)

                    [a,b,c,d,e,f,g,h]=[a*q+b,a,c*q+d,c,e*q+f,e,g*q+h,g]
            
                    y=1/(y-q)
                    iy+=1
    
    return exp


def JP_MCF(x,y,lim):
    [qa,qb]=[floor(x),floor(y)]
    
    [a,b]=[[qa],[qb]]

    for i in range(lim):
        [x,y]=[1/(y-qb),(x-qa)/(y-qb)]
        [qa,qb]=[floor(x),floor(y)]
        
        a.append(qa)
        b.append(qb)

    return [a,b]


#given the MCF of (x,y) compute the MCF of ((c11x+c12y+c13)/(c31x+c32y+c33),(c21x+c22y+c23)/(c31x+c32y+c33))
#where we declare the matrix of coefficients as C=matrix(RR,[[c11,c12,c13],[c21,c22,c23],[c31,c32,c33]])

def Jacobi_Gosper(C,x,y,lim):
    exp=[[],[]]

    
    for i in range(lim):
        print(i)
        if C[2][0]!=0 and C[2][1]!=0 and C[2][2]!=0:
            if floor(C[0][0]/C[2][0])==floor(C[0][1]/C[2][1])==floor(C[0][2]/C[2][2]) and floor(C[1][0]/C[2][0])==floor(C[1][1]/C[2][1])==floor(C[1][2]/C[2][2]) and sign(C[2][0])==sign(C[2][1])==sign(C[2][2]):
                qtrans1=floor(C[0][0]/C[2][0])
                qtrans2=floor(C[1][0]/C[2][0])
                exp[0].append(qtrans1)
                exp[1].append(qtrans2)

                trans_mat=matrix(RR,[[0,0,1],[1,0,-qtrans1],[0,1,-qtrans2]])
                C=trans_mat*C
                
            else:
                [q1,q2]=[floor(x),floor(y)]

                trans_mat=matrix(RR,[[q1,1,0],[q2,0,1],[1,0,0]])
                C=C*trans_mat

                [x,y]=[1/(y-q2),(x-q1)/(y-q2)]


        else:
            [q1,q2]=[floor(x),floor(y)]

            trans_mat=matrix(RR,[[q1,1,0],[q2,0,1],[1,0,0]])
            C=C*trans_mat

            [x,y]=[1/(y-q2),(x-q1)/(y-q2)]

    return exp


def seq_Jacobi_Gosper(C,a,b,lim):
    exp=[[],[]]

    j=0

    for i in range(lim):
        if C[2][0]!=0 and C[2][1]!=0 and C[2][2]!=0:
            if floor(C[0][0]/C[2][0])==floor(C[0][1]/C[2][1])==floor(C[0][2]/C[2][2]) and floor(C[1][0]/C[2][0])==floor(C[1][1]/C[2][1])==floor(C[1][2]/C[2][2]) and sign(C[2][0])==sign(C[2][1])==sign(C[2][2]):
                qtrans1=floor(C[0][0]/C[2][0])
                qtrans2=floor(C[1][0]/C[2][0])
                exp[0].append(qtrans1)
                exp[1].append(qtrans2)

                trans_mat=matrix(RR,[[0,0,1],[1,0,-qtrans1],[0,1,-qtrans2]])
                C=trans_mat*C
                
            else:
                [q1,q2]=[a[j],b[j]]

                trans_mat=matrix(RR,[[q1,1,0],[q2,0,1],[1,0,0]])
                C=C*trans_mat

                j+=1

        else:
            [q1,q2]=[a[j],b[j]]

            trans_mat=matrix(RR,[[q1,1,0],[q2,0,1],[1,0,0]])
            C=C*trans_mat

            j+=1

    return exp




#given the MCF of (x1,x1) and (y1,y2) compute the MCF of (C1(x1,x2)(y1,y2)/C3(x1,x2)(y1,y2),C2(x1,x2)(y1,y2)/C3(x1,x2)(y1,y2), where C1, C2,C3 are 3x3 matrices defining a bilinear form in the 9 variables x1y1,x1y2,x2y1,x2y2,x1,y1,x2,y2,1
#we declare the matrices as C=matrix(RR,[[c11,c12,c13],[c21,c22,c23],[c31,c32,c33]])

def old_Jacobi_Gosper2(C1,C2,C3,x1,x2,y1,y2,lim):
    exp=[[],[]]

    C1=matrix(RR,C1)
    C2=matrix(RR,C2)
    C3=matrix(RR,C3)
    [ix,iy]=[0,0]
    ninput=0

    for i in range(lim):
        if C1[2][0]!=0 and C1[2][1]!=0 and C1[2][2]!=0 and C2[2][0]!=0 and C2[2][1]!=0 and C2[2][2]!=0:
            if floor(C1[0][0]/C3[0][0])==floor(C1[0][1]/C3[0][1])==floor(C1[0][2]/C3[0][2])==floor(C1[1][0]/C3[1][0])==floor(C1[1][1]/C3[1][1])==floor(C1[1][2]/C3[1][2])==floor(C1[2][0]/C3[2][0])==floor(C1[2][1]/C3[2][1])==floor(C1[2][2]/C3[2][2]) and floor(C2[0][0]/C3[0][0])==floor(C2[0][1]/C3[0][1])==floor(C2[0][2]/C3[0][2])==floor(C2[1][0]/C3[1][0])==floor(C2[1][1]/C3[1][1])==floor(C2[1][2]/C3[1][2])==floor(C2[2][0]/C3[2][0])==floor(C2[2][1]/C3[2][1])==floor(C2[2][2]/C3[2][2]) and sign(C3[0][0])==sign(C3[0][1])==sign(C3[0][2])==sign(C3[1][0])==sign(C3[1][1])==sign(C3[1][2])==sign(C3[2][0])==sign(C3[2][1])==sign(C3[2][2]):
                [qtrans1,qtrans2]=[floor(C1[0][0]/C3[0][0]),floor(C2[0][0]/C3[0][0])]
                
                exp[0].append(qtrans1)
                exp[1].append(qtrans2)
               
                [C1,C2,C3]=[C3,C1-qtrans1*C3,C2-qtrans2*C3]
                
                
            else:
                ninput+=1
                print("The inputs are now:",ninput)

                if ix==iy:
                    [q1,q2]=[floor(x1),floor(x2)]
                    
                    trans_mat=matrix(RR,[[q1,q2,1],[1,0,0],[0,1,0]])
               
               
                    C1= trans_mat*C1
                   
                    
                    C2= trans_mat*C2
                    
                    C3= trans_mat*C3

                    [x1,x2]=[1/(x2-q2),(x1-q1)/(x2-q2)]

                    C1=C1.transpose()
                    C2=C2.transpose()
                    C3=C3.transpose()
                    
                    ix+=1

                else:
                    [q1,q2]=[floor(y1),floor(y2)]
                    
                    trans_mat=matrix(RR,[[q1,q2,1],[1,0,0],[0,1,0]])

                    C1= trans_mat*C1

                    C2= trans_mat*C2
                    
                    C3= trans_mat*C3
            
                    [y1,y2]=[1/(y2-q2),(y1-q1)/(y2-q2)]

                    C1=C1.transpose()
                    C2=C2.transpose()
                    C3=C3.transpose()

                    iy+=1

        else:
            ninput+=1
            print("The inputs are now:",ninput)
            if ix==iy:
                [q1,q2]=[floor(x1),floor(x2)]
                
                trans_mat=matrix(RR,[[q1,q2,1],[1,0,0],[0,1,0]])
                C1= trans_mat*C1

                C2= trans_mat*C2
                
                C3= trans_mat*C3

                [x1,x2]=[1/(x2-q2),(x1-q1)/(x2-q2)]

                C1=C1.transpose()
                C2=C2.transpose()
                C3=C3.transpose()
                
                ix+=1

            else:
                [q1,q2]=[floor(y1),floor(y2)]
                
                trans_mat=matrix(RR,[[q1,q2,1],[1,0,0],[0,1,0]])

                C1= trans_mat*C1

                C2= trans_mat*C2
                
                C3= trans_mat*C3

                C1=C1.transpose()
                C2=C2.transpose()
                C3=C3.transpose()
        
                [y1,y2]=[1/(y2-q2),(y1-q1)/(y2-q2)]

                iy+=1
                
    return exp

def Jacobi_Gosper2(C1,C2,C3,x1,x2,y1,y2,lim):
    exp=[[],[]]

    C1=matrix(RR,C1)
    C2=matrix(RR,C2)
    C3=matrix(RR,C3)
    [ix,iy]=[0,0]
    ninput=0

    for i in range(lim):
        if C3[0][0]!=0 and C3[0][1]!=0 and C3[0][2]!=0 and C3[1][0]!=0 and C3[1][1]!=0 and C3[1][2]!=0 and C3[2][0]!=0 and C3[2][1]!=0 and C3[2][2]!=0:
            if floor(C1[0][0]/C3[0][0])==floor(C1[0][1]/C3[0][1])==floor(C1[0][2]/C3[0][2])==floor(C1[1][0]/C3[1][0])==floor(C1[1][1]/C3[1][1])==floor(C1[1][2]/C3[1][2])==floor(C1[2][0]/C3[2][0])==floor(C1[2][1]/C3[2][1])==floor(C1[2][2]/C3[2][2]) and floor(C2[0][0]/C3[0][0])==floor(C2[0][1]/C3[0][1])==floor(C2[0][2]/C3[0][2])==floor(C2[1][0]/C3[1][0])==floor(C2[1][1]/C3[1][1])==floor(C2[1][2]/C3[1][2])==floor(C2[2][0]/C3[2][0])==floor(C2[2][1]/C3[2][1])==floor(C2[2][2]/C3[2][2]) and sign(C3[0][0])==sign(C3[0][1])==sign(C3[0][2])==sign(C3[1][0])==sign(C3[1][1])==sign(C3[1][2])==sign(C3[2][0])==sign(C3[2][1])==sign(C3[2][2]):
                [qtrans1,qtrans2]=[floor(C1[0][0]/C3[0][0]),floor(C2[0][0]/C3[0][0])]
                
                exp[0].append(qtrans1)
                exp[1].append(qtrans2)
               
                [C1,C2,C3]=[C3,C1-qtrans1*C3,C2-qtrans2*C3]
                
                
            else:
                ninput+=1
                print("The inputs are now:",ninput)

                if ix==iy:
                    [q1,q2]=[floor(x1),floor(x2)]
                    
                    trans_mat=matrix(RR,[[q1,q2,1],[1,0,0],[0,1,0]])
               
               
                    C1= trans_mat*C1
                   
                    
                    C2= trans_mat*C2
                    
                    C3= trans_mat*C3

                    [x1,x2]=[1/(x2-q2),(x1-q1)/(x2-q2)]

                    C1=C1.transpose()
                    C2=C2.transpose()
                    C3=C3.transpose()
                    
                    ix+=1

                else:
                    [q1,q2]=[floor(y1),floor(y2)]
                    
                    trans_mat=matrix(RR,[[q1,q2,1],[1,0,0],[0,1,0]])

                    C1= trans_mat*C1

                    C2= trans_mat*C2
                    
                    C3= trans_mat*C3
            
                    [y1,y2]=[1/(y2-q2),(y1-q1)/(y2-q2)]

                    C1=C1.transpose()
                    C2=C2.transpose()
                    C3=C3.transpose()

                    iy+=1

        else:
            ninput+=1
            print("The inputs are now:",ninput)
            if ix==iy:
                [q1,q2]=[floor(x1),floor(x2)]
                
                trans_mat=matrix(RR,[[q1,q2,1],[1,0,0],[0,1,0]])
                C1= trans_mat*C1

                C2= trans_mat*C2
                
                C3= trans_mat*C3

                [x1,x2]=[1/(x2-q2),(x1-q1)/(x2-q2)]

                C1=C1.transpose()
                C2=C2.transpose()
                C3=C3.transpose()
                
                ix+=1

            else:
                [q1,q2]=[floor(y1),floor(y2)]
                
                trans_mat=matrix(RR,[[q1,q2,1],[1,0,0],[0,1,0]])

                C1= trans_mat*C1

                C2= trans_mat*C2
                
                C3= trans_mat*C3

                C1=C1.transpose()
                C2=C2.transpose()
                C3=C3.transpose()
        
                [y1,y2]=[1/(y2-q2),(y1-q1)/(y2-q2)]

                iy+=1
                
    return exp


def seq_Jacobi_Gosper2(C1,C2,C3,x1,x2,y1,y2,lim):
    exp=[[],[]]

    C1=matrix(RR,C1)
    C2=matrix(RR,C2)
    C3=matrix(RR,C3)
    [ix,iy]=[0,0]
    ninput=0

    for i in range(lim):
        print(C1)
        print('SPACE')
        print(C2)
        print('SPACE')
        print(C3)
        print('SPACEfine')

        if C3[0][0]!=0 and C3[0][1]!=0 and C3[0][2]!=0 and C3[1][0]!=0 and C3[1][1]!=0 and C3[1][2]!=0 and C3[2][0]!=0 and C3[2][1]!=0 and C3[2][2]!=0:
            if floor(C1[0][0]/C3[0][0])==floor(C1[0][1]/C3[0][1])==floor(C1[0][2]/C3[0][2])==floor(C1[1][0]/C3[1][0])==floor(C1[1][1]/C3[1][1])==floor(C1[1][2]/C3[1][2])==floor(C1[2][0]/C3[2][0])==floor(C1[2][1]/C3[2][1])==floor(C1[2][2]/C3[2][2]) and floor(C2[0][0]/C3[0][0])==floor(C2[0][1]/C3[0][1])==floor(C2[0][2]/C3[0][2])==floor(C2[1][0]/C3[1][0])==floor(C2[1][1]/C3[1][1])==floor(C2[1][2]/C3[1][2])==floor(C2[2][0]/C3[2][0])==floor(C2[2][1]/C3[2][1])==floor(C2[2][2]/C3[2][2]) and sign(C3[0][0])==sign(C3[0][1])==sign(C3[0][2])==sign(C3[1][0])==sign(C3[1][1])==sign(C3[1][2])==sign(C3[2][0])==sign(C3[2][1])==sign(C3[2][2]):
                [qtrans1,qtrans2]=[floor(C1[0][0]/C3[0][0]),floor(C2[0][0]/C3[0][0])]
                

                print('IM DOING OUTPUT!')
                exp[0].append(qtrans1)
                exp[1].append(qtrans2)
               
                [C1,C2,C3]=[C3,C1-qtrans1*C3,C2-qtrans2*C3]
                
                
            else:
                ninput+=1
                print("The inputs are now:",ninput)

                if ix==iy:
                    [q1,q2]=[x1[ix],x2[ix]]
                    
                    trans_mat=matrix(RR,[[q1,q2,1],[1,0,0],[0,1,0]])

                    C1= trans_mat*C1
                    C2= trans_mat*C2
                    C3= trans_mat*C3

                    C1=C1.transpose()
                    C2=C2.transpose()
                    C3=C3.transpose()
                    
                    ix+=1


                else:
                    [q1,q2]=[y1[iy],y2[iy]]
                    
                    trans_mat=matrix(RR,[[q1,q2,1],[1,0,0],[0,1,0]])
                    
                    
                    C1= trans_mat*C1

                    C2= trans_mat*C2
                    
                    C3= trans_mat*C3
            

                    C1=C1.transpose()
                    C2=C2.transpose()
                    C3=C3.transpose()

                    iy+=1
                    


        else:
            ninput+=1
            print("The inputs are now:",ninput)

            if ix==iy:
                [q1,q2]=[x1[ix],x2[ix]]
                
                trans_mat=matrix(RR,[[q1,q2,1],[1,0,0],[0,1,0]])
            
            
                C1= trans_mat*C1
                C2= trans_mat*C2
                C3= trans_mat*C3

                C1=C1.transpose()
                C2=C2.transpose()
                C3=C3.transpose()
                
                ix+=1


            else:
                [q1,q2]=[y1[iy],y2[iy]]
                
                trans_mat=matrix(RR,[[q1,q2,1],[1,0,0],[0,1,0]])

                C1= trans_mat*C1

                C2= trans_mat*C2
                
                C3= trans_mat*C3
        

                C1=C1.transpose()
                C2=C2.transpose()
                C3=C3.transpose()

                iy+=1
                
    return exp



#given d>0 non-cube, from the continued fraction (d^(1/3),d^(2/3)), we expand (..,...) and collect data of:
#number of inputs before outputs at each step, plot the growing graph of inputs  (n on the x-axis, inputs on the y-axis), it is an increasing function


def count_inputs(d,lim):
    exp=[[],[]]

    [x,y]=[d^(1/3),d^(2/3)]
    C=matrix(RR,[[1,0,0],[0,1,0],[0,0,1]])

    inputs_array=[]
    total_inputs_array=[]
    total_inputs=0
    total_outputs=0
    input_counter=0

    
    for i in range(lim):
        print("Step:",i)
        if C[2][0]!=0 and C[2][1]!=0 and C[2][2]!=0:
            if floor(C[0][0]/C[2][0])==floor(C[0][1]/C[2][1])==floor(C[0][2]/C[2][2]) and floor(C[1][0]/C[2][0])==floor(C[1][1]/C[2][1])==floor(C[1][2]/C[2][2]) and sign(C[2][0])==sign(C[2][1])==sign(C[2][2]):
                qtrans1=floor(C[0][0]/C[2][0])
                qtrans2=floor(C[1][0]/C[2][0])
                exp[0].append(qtrans1)
                exp[1].append(qtrans2)

                total_outputs+=1                           #total number of outputs
                inputs_array.append(input_counter)         #inputs required for this step
                total_inputs_array.append(total_inputs)    #total inputs required until this output
                input_counter=0                            

                trans_mat=matrix(RR,[[0,0,1],[1,0,-qtrans1],[0,1,-qtrans2]])
                C=trans_mat*C
                
            else:
                [q1,q2]=[floor(x),floor(y)]

                input_counter+=1
                total_inputs+=1

                trans_mat=matrix(RR,[[q1,1,0],[q2,0,1],[1,0,0]])
                C=C*trans_mat

                [x,y]=[1/(y-q2),(x-q1)/(y-q2)]

        else:
            [q1,q2]=[floor(x),floor(y)]

            input_counter+=1
            total_inputs+=1

            trans_mat=matrix(RR,[[q1,1,0],[q2,0,1],[1,0,0]])
            C=C*trans_mat
            
            [x,y]=[1/(y-q2),(x-q1)/(y-q2)]

    print("Total inputs:",total_inputs)
    print("Total inputs array:",total_inputs_array)
    print("Inputs array:",inputs_array)
    print("Total outputs:",total_outputs)

    return total_inputs_array


def seq_count_inputs(C,a,b,lim):
    exp=[[],[]]

    
    C=matrix(RR,[[3,0,0],[0,-2,0],[0,0,6]])

    print(C)

    u=0

    inputs_array=[]
    total_inputs_array=[]
    total_inputs=0
    total_outputs=0
    input_counter=0
    
    for i in range(lim):
        print("Step:",i)
        if C[2][0]!=0 and C[2][1]!=0 and C[2][2]!=0:
            if floor(C[0][0]/C[2][0])==floor(C[0][1]/C[2][1])==floor(C[0][2]/C[2][2]) and floor(C[1][0]/C[2][0])==floor(C[1][1]/C[2][1])==floor(C[1][2]/C[2][2]) and sign(C[2][0])==sign(C[2][1])==sign(C[2][2]):
                qtrans1=floor(C[0][0]/C[2][0])
                qtrans2=floor(C[1][0]/C[2][0])
                exp[0].append(qtrans1)
                exp[1].append(qtrans2)

                total_outputs+=1                           #total number of outputs
                inputs_array.append(input_counter)         #inputs required for this step
                total_inputs_array.append(total_inputs)    #total inputs required until this output
                input_counter=0                            

                trans_mat=matrix(RR,[[0,0,1],[1,0,-qtrans1],[0,1,-qtrans2]])
                C=trans_mat*C
                
            else:
                [q1,q2]=[a[u],b[u]]

                trans_mat=matrix(RR,[[q1,1,0],[q2,0,1],[1,0,0]])
                C=C*trans_mat

                u+=1

                input_counter+=1
                total_inputs+=1

        else:
            [q1,q2]=[a[u],b[u]]

            trans_mat=matrix(RR,[[q1,1,0],[q2,0,1],[1,0,0]])
            C=C*trans_mat

            u+=1

            input_counter+=1
            total_inputs+=1

    print("Total inputs:",total_inputs)
    print("Total inputs array:",total_inputs_array)
    print("Inputs array:",inputs_array)
    print("Total outputs:",total_outputs)

    return total_inputs_array


#the following function, fixing the number of outputs, gives back the total number of inputs requiredß

def count_inputs_upto(d,out_lim):
    exp=[[],[]]

    [x,y]=[d^(1/3),d^(2/3)]
    C=matrix(RR,[[3,5,0],[5,3,0],[1,0,2]])
    print(C)
    inputs_array=[]
    total_inputs_array=[]
    total_inputs=0
    total_outputs=0
    input_counter=0

    i=0
    
    while total_outputs<out_lim:
        print("Step:",i)
        if C[2][0]!=0 and C[2][1]!=0 and C[2][2]!=0:
            if floor(C[0][0]/C[2][0])==floor(C[0][1]/C[2][1])==floor(C[0][2]/C[2][2]) and floor(C[1][0]/C[2][0])==floor(C[1][1]/C[2][1])==floor(C[1][2]/C[2][2]) and sign(C[2][0])==sign(C[2][1])==sign(C[2][2]):
                qtrans1=floor(C[0][0]/C[2][0])
                qtrans2=floor(C[1][0]/C[2][0])
                exp[0].append(qtrans1)
                exp[1].append(qtrans2)

                total_outputs+=1                           #total number of outputs
                inputs_array.append(input_counter)         #inputs required for this step
                total_inputs_array.append(total_inputs)    #total inputs required until this output
                input_counter=0                            

                trans_mat=matrix(RR,[[0,0,1],[1,0,-qtrans1],[0,1,-qtrans2]])
                C=trans_mat*C
                
            else:
                [q1,q2]=[floor(x),floor(y)]

                input_counter+=1
                total_inputs+=1

                trans_mat=matrix(RR,[[q1,1,0],[q2,0,1],[1,0,0]])
                C=C*trans_mat

                [x,y]=[1/(y-q2),(x-q1)/(y-q2)]

        else:
            [q1,q2]=[floor(x),floor(y)]

            input_counter+=1
            total_inputs+=1

            trans_mat=matrix(RR,[[q1,1,0],[q2,0,1],[1,0,0]])
            C=C*trans_mat
            
            [x,y]=[1/(y-q2),(x-q1)/(y-q2)]
        i+=1

    print("Total inputs:",total_inputs)
    print("Total inputs array:",total_inputs_array)
    print("Inputs array:",inputs_array)
    print("Total outputs:",total_outputs)

    return total_inputs_array


#for all non-cube integers up to D, I print the mean inputs required at each step to get an output, until out_lim outputs

def plot_mean_inputs_upto(D,out_lim):
    input_mean=[]
    inputs_sequences=[]
    for d in range(2,D):
        if (d^(1/3)).is_integer()==False:
            input_list=count_inputs_upto(d,out_lim)
            inputs_sequences.append(input_list)
        

    for k in range(out_lim):
        mean=0
        for i in range(len(inputs_sequences)):
            mean+=inputs_sequences[i][k]
        ll=len(inputs_sequences)
        mean=mean/ll
        input_mean.append(RR(mean))

    print("The sequences:",inputs_sequences)

    print("The mean:")
    return input_mean

def plot_inputs(d,lim):

    C=matrix(RR,[[2,0,0],[0,2,0],[0,0,1]])
    inp_ar=count_inputs(d,lim)
    len1=len(inp_ar)

    coord_per1= [(i,inp_ar[i]) for i in range(len1)]
    G1=list_plot(coord_per1, gridlines=True, plotjoined=True, marker='.',color='blue',thickness=2,legend_label='Inputs')

    Y=G1+text("Outputs",(10,-2),color='black',fontsize=14)+text("Inputs required",(-2,10),color='black',fontsize=14, rotation="vertical")
    return Y


#SINGLE m=1

def uni_count_random_inputs_upto(C,out_lim,pq_lim):

    exp=[]

    starting_cf=[[],[]]

    print(C)

    [a,b,c,d]=[C[0][0],C[0][1],C[1][0],C[1][1]]

    print(a,b,c,d)

    i=0

    inputs_array=[]
    total_inputs_array=[]
    total_inputs=0
    total_outputs=0
    input_counter=0

    while total_outputs<out_lim:
        print("Step:",i)
        if c!=0 and (c+d)!=0:
            if floor(a/c)==floor((a+b)/(c+d)) and sign(c)==sign(c+d):
                qtrans=floor(a/c)
                exp.append(qtrans)

                total_outputs+=1                           #total number of outputs
                inputs_array.append(input_counter)         #inputs required for this step
                total_inputs_array.append(total_inputs)    #total inputs required until this output
                input_counter=0                            

                [a,b,c,d]=[c,d,a-c*qtrans,b-d*qtrans]
                
            else:
                q=randint(1,pq_lim)

                [a,b,c,d]=[a*q+b,a,c*q+d,c]
                
                input_counter+=1
                total_inputs+=1

        else:
                q=randint(1,pq_lim)

                [a,b,c,d]=[a*q+b,a,c*q+d,c]
                
                input_counter+=1
                total_inputs+=1

        i+=1

    print("Total inputs:",total_inputs)
    print("Total inputs array:",total_inputs_array)
    print("Inputs array:",inputs_array)
    print("Total outputs:",total_outputs)

    print("Initial expansion is:",starting_cf)
    print("Final expansion:",exp)

    return total_inputs_array


def uni_plot_random_inputs_without_mean_upto(out_lim,pq_lim,coeff_lim,N,col):

    inputs_ar=[]

    c1=random_vector(2,coeff_lim)
    c2=random_vector(2,coeff_lim)

    C=matrix(RR,[c1,c2])

    while C.determinant()==0:
        c1=random_vector(2,coeff_lim)
        c2=random_vector(2,coeff_lim)

        C=matrix(RR,[c1,c2])


    for i in range(N):
        total_inp_ar=uni_count_random_inputs_upto(C,out_lim,pq_lim)
        inputs_ar.append(total_inp_ar)

    
    seq_to_plot=[]

    for i in range(len(inputs_ar)):
        plot_to_app=[]
        for j in range(out_lim):
            plot_to_app.append([j+1,inputs_ar[i][j]])
        seq_to_plot.append(plot_to_app)

    G=list_plot(seq_to_plot[0], gridlines=True, plotjoined=True,color=col,thickness=1)

    for i in range(1,len(inputs_ar)):
        G=G+list_plot(seq_to_plot[i], gridlines=True, plotjoined=True,color=col,thickness=1)
    
    return G


#BILINEAR m=1

    exp=[]

    [ix,iy]=[0,0]

    for i in range(lim):
        if e!=0 and (e+f)!=0 and (e+g)!=0 and (e+f+g+h)!=0:
            if floor(a/e)==floor((a+b)/(e+f))==floor((a+c)/(e+g))==floor((a+b+c+d)/(e+f+g+h)) and  sign(e)==sign(e+f)==sign(e+g):
                qtrans=floor(a/e)
                exp.append(qtrans)

                [a,b,c,d,e,f,g,h]=[e,f,g,h,a-e*qtrans,b-f*qtrans,c-g*qtrans,d-h*qtrans]
                
            else:
                if ix==iy:
                    q=floor(x)

                    [a,b,c,d,e,f,g,h]=[a*q+c,b*q+d,a,b,e*q+g,f*q+h,e,f]

                    x=1/(x-q)
                    ix+=1

                else:
                    q=floor(y)

                    [a,b,c,d,e,f,g,h]=[a*q+b,a,c*q+d,c,e*q+f,e,g*q+h,g]
            
                    y=1/(y-q)
                    iy+=1

        else:
                if ix==iy:
                    q=floor(x)

                    [a,b,c,d,e,f,g,h]=[a*q+c,b*q+d,a,b,e*q+g,f*q+h,e,f]

                    x=1/(x-q)
                    ix+=1

                else:
                    q=floor(y)

                    [a,b,c,d,e,f,g,h]=[a*q+b,a,c*q+d,c,e*q+f,e,g*q+h,g]
            
                    y=1/(y-q)
                    iy+=1
    
    return exp

def uni_biquad_count_random_inputs_upto(C1,C2,out_lim,pq_lim):

    exp=[]

    [a,b,c,d,e,f,g,h]=[C1[0][0],C1[0][1],C1[1][0],C1[1][1],C2[0][0],C2[0][1],C2[1][0],C2[1][1]]

    [ix,iy]=[0,0]

    x1=[]
    y1=[]

    for i in range(10*out_lim):
        com1=randint(1,pq_lim)
        x1.append(com1)

    for i in range(10*out_lim):
        com1=randint(1,pq_lim)
        y1.append(com1)
    i=0

    inputs_array=[]
    total_inputs_array=[]
    total_inputs=0
    total_outputs=0
    input_counter=0

    while total_outputs<out_lim:
        print("Step:",i)

        if e!=0 and (e+f)!=0 and (e+g)!=0 and (e+f+g+h)!=0:
            if floor(a/e)==floor((a+b)/(e+f))==floor((a+c)/(e+g))==floor((a+b+c+d)/(e+f+g+h)) and  sign(e)==sign(e+f)==sign(e+g):
                qtrans=floor(a/e)
                exp.append(qtrans)

                [a,b,c,d,e,f,g,h]=[e,f,g,h,a-e*qtrans,b-f*qtrans,c-g*qtrans,d-h*qtrans]
                
                exp.append(qtrans)

                total_outputs+=1                           #total number of outputs
                inputs_array.append(input_counter)         #inputs required for this step
                total_inputs_array.append(total_inputs)    #total inputs required until this output
                input_counter=0                            


            else:
                if ix==iy:
                    q=x1[ix]

                    [a,b,c,d,e,f,g,h]=[a*q+c,b*q+d,a,b,e*q+g,f*q+h,e,f]

                    input_counter+=1
                    total_inputs+=1

                    ix+=1

                else:
                    q=y1[iy]

                    [a,b,c,d,e,f,g,h]=[a*q+b,a,c*q+d,c,e*q+f,e,g*q+h,g]

                    input_counter+=1
                    total_inputs+=1
            
               
                    iy+=1
                    


        else:
            if ix==iy:
                    q=x1[ix]

                    [a,b,c,d,e,f,g,h]=[a*q+c,b*q+d,a,b,e*q+g,f*q+h,e,f]

                    input_counter+=1
                    total_inputs+=1

                    ix+=1

            else:
                q=y1[iy]

                [a,b,c,d,e,f,g,h]=[a*q+b,a,c*q+d,c,e*q+f,e,g*q+h,g]

                input_counter+=1
                total_inputs+=1
        
            
                iy+=1
            

    print("Total inputs:",total_inputs)
    print("Total inputs array:",total_inputs_array)
    print("Inputs array:",inputs_array)
    print("Total outputs:",total_outputs)

    return total_inputs_array


def uni_biquad_plot_random_inputs_without_mean_upto(out_lim,pq_lim,coeff_lim,N,col):

    inputs_ar=[]

    c11=random_vector(2,coeff_lim)
    c12=random_vector(2,coeff_lim)

    c21=random_vector(2,coeff_lim)
    c22=random_vector(2,coeff_lim)

    C1=matrix(RR,[c11,c12])
    C2=matrix(RR,[c21,c22])
    

    while C1.determinant()==0:
        c11=random_vector(2,coeff_lim)
        c12=random_vector(2,coeff_lim)

        C1=matrix(RR,[c11,c12])

    while C2.determinant()==0:
        c21=random_vector(2,coeff_lim)
        c22=random_vector(2,coeff_lim)

        C2=matrix(RR,[c21,c22])


    for i in range(N):
        total_inp_ar=uni_biquad_count_random_inputs_upto(C1,C2,out_lim,pq_lim)
        inputs_ar.append(total_inp_ar)

    
    seq_to_plot=[]

    for i in range(len(inputs_ar)):
        plot_to_app=[]
        for j in range(out_lim):
            plot_to_app.append([j+1,inputs_ar[i][j]])
        seq_to_plot.append(plot_to_app)

    G=list_plot(seq_to_plot[0], gridlines=True, plotjoined=True,color=col,thickness=1)

    for i in range(1,len(inputs_ar)):
        G=G+list_plot(seq_to_plot[i], gridlines=True, plotjoined=True,color=col,thickness=1)
    
    return G





#SINGLE m=2

def count_random_inputs_upto(C,out_lim,pq_lim):

    exp=[[],[]]

    starting_cf=[[],[]]

    print(C)

    i=0

    inputs_array=[]
    total_inputs_array=[]
    total_inputs=0
    total_outputs=0
    input_counter=0

    while total_outputs<out_lim:
        print("Step:",i)
        if C[2][0]!=0 and C[2][1]!=0 and C[2][2]!=0:
            if floor(C[0][0]/C[2][0])==floor(C[0][1]/C[2][1])==floor(C[0][2]/C[2][2]) and floor(C[1][0]/C[2][0])==floor(C[1][1]/C[2][1])==floor(C[1][2]/C[2][2]) and sign(C[2][0])==sign(C[2][1])==sign(C[2][2]):
                qtrans1=floor(C[0][0]/C[2][0])
                qtrans2=floor(C[1][0]/C[2][0])
                exp[0].append(qtrans1)
                exp[1].append(qtrans2)

                total_outputs+=1                           #total number of outputs
                inputs_array.append(input_counter)         #inputs required for this step
                total_inputs_array.append(total_inputs)    #total inputs required until this output
                input_counter=0                            

                trans_mat=matrix(RR,[[0,0,1],[1,0,-qtrans1],[0,1,-qtrans2]])
                C=trans_mat*C
                
            else:

                q1=randint(1,pq_lim)
                q2=randint(0,q1)
                if i==0:
                    q2=randint(0,pq_lim)
                
                starting_cf[0].append(q1)
                starting_cf[1].append(q2)
  
                input_counter+=1
                total_inputs+=1

                trans_mat=matrix(RR,[[q1,1,0],[q2,0,1],[1,0,0]])
                C=C*trans_mat

        else:
            q1=randint(1,pq_lim)
            q2=randint(0,q1)
            if i==0:
                q2=randint(0,pq_lim)

            starting_cf[0].append(q1)
            starting_cf[1].append(q2)

            input_counter+=1
            total_inputs+=1

            trans_mat=matrix(RR,[[q1,1,0],[q2,0,1],[1,0,0]])
            C=C*trans_mat

        i+=1

    print("Total inputs:",total_inputs)
    print("Total inputs array:",total_inputs_array)
    print("Inputs array:",inputs_array)
    print("Total outputs:",total_outputs)

    print("Initial expansion is:",starting_cf)
    print("Final expansion:",exp)

    return total_inputs_array


def plot_random_inputs_upto(out_lim,pq_lim,coeff_lim,N,col):

    inputs_ar=[]

    for i in range(N):
        c1=random_vector(3,coeff_lim)
        c2=random_vector(3,coeff_lim)
        c3=random_vector(3,coeff_lim)

        C=matrix(RR,[c1,c2,c3])

        while C.determinant()==0:
            c1=random_vector(3,coeff_lim)
            c2=random_vector(3,coeff_lim)
            c3=random_vector(3,coeff_lim)

            C=matrix(RR,[c1,c2,c3])

        total_inp_ar=count_random_inputs_upto(C,out_lim,pq_lim)
        inputs_ar.append(total_inp_ar)

    mean_inputs_ar=[]

    for k in range(out_lim):
        mean=0
        for i in range(N):
            mean+=inputs_ar[i][k]
        
        mean=mean/N
        mean_inputs_ar.append(RR(mean))

    points_for_plot=[]
    for i in range(out_lim):
        points_for_plot.append([i+1,mean_inputs_ar[i]])
    print(out_lim)
    G=list_plot(points_for_plot, gridlines=True, plotjoined=True,color=col,thickness=1)
   
    return G


def plot_random_inputs_without_mean_upto(out_lim,pq_lim,coeff_lim,N,col):

    inputs_ar=[]

    c1=random_vector(3,coeff_lim)
    c2=random_vector(3,coeff_lim)
    c3=random_vector(3,coeff_lim)

    C=matrix(RR,[c1,c2,c3])

    while C.determinant()==0:
        c1=random_vector(3,coeff_lim)
        c2=random_vector(3,coeff_lim)
        c3=random_vector(3,coeff_lim)

        C=matrix(RR,[c1,c2,c3])


    for i in range(N):
        total_inp_ar=count_random_inputs_upto(C,out_lim,pq_lim)
        inputs_ar.append(total_inp_ar)

    
    seq_to_plot=[]

    for i in range(len(inputs_ar)):
        plot_to_app=[]
        for j in range(out_lim):
            plot_to_app.append([j+1,inputs_ar[i][j]])
        seq_to_plot.append(plot_to_app)

    G=list_plot(seq_to_plot[0], gridlines=True, plotjoined=True,color=col,thickness=1)

    for i in range(1,len(inputs_ar)):
        G=G+list_plot(seq_to_plot[i], gridlines=True, plotjoined=True,color=col,thickness=1)
    
    return G



#BILINEAR m=2


def biquad_count_random_inputs_upto(C1,C2,C3,out_lim,pq_lim):

    exp=[[],[]]

    starting_cf=[[],[]]

    [ix,iy]=[0,0]
    ninput=0

    x1=[]
    x2=[]
    y1=[]
    y2=[]

    for i in range(10*out_lim):
        com1=randint(1,pq_lim)
        com2=randint(0,com1)
        x1.append(com1)
        x2.append(com2)

    for i in range(10*out_lim):
        com1=randint(1,pq_lim)
        com2=randint(0,com1)
        y1.append(com1)
        y2.append(com2)

    i=0

    inputs_array=[]
    total_inputs_array=[]
    total_inputs=0
    total_outputs=0
    input_counter=0

    while total_outputs<out_lim:
        print("Step:",i)

        if C3[0][0]!=0 and C3[0][1]!=0 and C3[0][2]!=0 and C3[1][0]!=0 and C3[1][1]!=0 and C3[1][2]!=0 and C3[2][0]!=0 and C3[2][1]!=0 and C3[2][2]!=0:
            if floor(C1[0][0]/C3[0][0])==floor(C1[0][1]/C3[0][1])==floor(C1[0][2]/C3[0][2])==floor(C1[1][0]/C3[1][0])==floor(C1[1][1]/C3[1][1])==floor(C1[1][2]/C3[1][2])==floor(C1[2][0]/C3[2][0])==floor(C1[2][1]/C3[2][1])==floor(C1[2][2]/C3[2][2]) and floor(C2[0][0]/C3[0][0])==floor(C2[0][1]/C3[0][1])==floor(C2[0][2]/C3[0][2])==floor(C2[1][0]/C3[1][0])==floor(C2[1][1]/C3[1][1])==floor(C2[1][2]/C3[1][2])==floor(C2[2][0]/C3[2][0])==floor(C2[2][1]/C3[2][1])==floor(C2[2][2]/C3[2][2]) and sign(C3[0][0])==sign(C3[0][1])==sign(C3[0][2])==sign(C3[1][0])==sign(C3[1][1])==sign(C3[1][2])==sign(C3[2][0])==sign(C3[2][1])==sign(C3[2][2]):
                [qtrans1,qtrans2]=[floor(C1[0][0]/C3[0][0]),floor(C2[0][0]/C3[0][0])]
                
                exp[0].append(qtrans1)
                exp[1].append(qtrans2)

                total_outputs+=1                           #total number of outputs
                inputs_array.append(input_counter)         #inputs required for this step
                total_inputs_array.append(total_inputs)    #total inputs required until this output
                input_counter=0                            


               
                [C1,C2,C3]=[C3,C1-qtrans1*C3,C2-qtrans2*C3]

                
                
                
            else:
                ninput+=1
                print("The inputs are now:",ninput)

                if ix==iy:
                    [q1,q2]=[x1[ix],x2[ix]]
                    
                    trans_mat=matrix(RR,[[q1,q2,1],[1,0,0],[0,1,0]])
               
               
                    C1= trans_mat*C1
                    
                    C2= trans_mat*C2
                    
                    C3= trans_mat*C3

                    C1=C1.transpose()
                    C2=C2.transpose()
                    C3=C3.transpose()

                    input_counter+=1
                    total_inputs+=1

                    

                    
                    ix+=1

                else:
                    [q1,q2]=[y1[iy],y2[iy]]
                    
                    trans_mat=matrix(RR,[[q1,q2,1],[1,0,0],[0,1,0]])

                    C1= trans_mat*C1

                    C2= trans_mat*C2
                    
                    C3= trans_mat*C3

                    input_counter+=1
                    total_inputs+=1

            

                    C1=C1.transpose()
                    C2=C2.transpose()
                    C3=C3.transpose()

                    iy+=1
                    


        else:
            ninput+=1
            print("The inputs are now:",ninput)

            if ix==iy:
                [q1,q2]=[x1[ix],x2[ix]]
                
                trans_mat=matrix(RR,[[q1,q2,1],[1,0,0],[0,1,0]])
            
            
                C1= trans_mat*C1
                
                
                C2= trans_mat*C2
                
                C3= trans_mat*C3

                C1=C1.transpose()
                C2=C2.transpose()
                C3=C3.transpose()

                input_counter+=1
                total_inputs+=1

                
                ix+=1


            else:
                [q1,q2]=[y1[iy],y2[iy]]
                
                trans_mat=matrix(RR,[[q1,q2,1],[1,0,0],[0,1,0]])

                C1= trans_mat*C1

                C2= trans_mat*C2
                
                C3= trans_mat*C3

                input_counter+=1
                total_inputs+=1

        

                C1=C1.transpose()
                C2=C2.transpose()
                C3=C3.transpose()

                iy+=1
                

    print("Total inputs:",total_inputs)
    print("Total inputs array:",total_inputs_array)
    print("Inputs array:",inputs_array)
    print("Total outputs:",total_outputs)

    print("Initial expansion is:",starting_cf)
    print("Final expansion:",exp)

    return total_inputs_array


def biquad_plot_random_inputs_without_mean_upto(out_lim,pq_lim,coeff_lim,N,col):

    inputs_ar=[]

    c11=random_vector(3,coeff_lim)
    c12=random_vector(3,coeff_lim)
    c13=random_vector(3,coeff_lim)

    c21=random_vector(3,coeff_lim)
    c22=random_vector(3,coeff_lim)
    c23=random_vector(3,coeff_lim)

    c31=random_vector(3,coeff_lim)
    c32=random_vector(3,coeff_lim)
    c33=random_vector(3,coeff_lim)

    C1=matrix(RR,[c11,c12,c13])
    C2=matrix(RR,[c21,c22,c23])
    C3=matrix(RR,[c31,c32,c33])
    

    while C1.determinant()==0:
        c11=random_vector(3,coeff_lim)
        c12=random_vector(3,coeff_lim)
        c13=random_vector(3,coeff_lim)

        C1=matrix(RR,[c11,c12,c13])

    while C2.determinant()==0:
        c21=random_vector(3,coeff_lim)
        c22=random_vector(3,coeff_lim)
        c23=random_vector(3,coeff_lim)

        C2=matrix(RR,[c21,c22,c23])

    while C3.determinant()==0:
        c31=random_vector(3,coeff_lim)
        c32=random_vector(3,coeff_lim)
        c33=random_vector(3,coeff_lim)

        C3=matrix(RR,[c31,c32,c33])


    for i in range(N):
        total_inp_ar=biquad_count_random_inputs_upto(C1,C2,C3,out_lim,pq_lim)
        inputs_ar.append(total_inp_ar)

    
    seq_to_plot=[]

    for i in range(len(inputs_ar)):
        plot_to_app=[]
        for j in range(out_lim):
            plot_to_app.append([j+1,inputs_ar[i][j]])
        seq_to_plot.append(plot_to_app)

    G=list_plot(seq_to_plot[0], gridlines=True, plotjoined=True,color=col,thickness=1)

    for i in range(1,len(inputs_ar)):
        G=G+list_plot(seq_to_plot[i], gridlines=True, plotjoined=True,color=col,thickness=1)
    
    return G


#SINGLE m=3

#In this case, C is a 4x4 matrix

def tri_count_random_inputs_upto(C,out_lim,pq_lim):

    exp=[[],[],[]]

    starting_cf=[[],[],[]]

    print(C)

    i=0

    inputs_array=[]
    total_inputs_array=[]
    total_inputs=0
    total_outputs=0
    input_counter=0

    while total_outputs<out_lim:
        print("Step:",i)
        if C[3][0]!=0 and C[3][1]!=0 and C[3][2]!=0 and C[3][3]!=0:
            if floor(C[0][0]/C[3][0])==floor(C[0][1]/C[3][1])==floor(C[0][2]/C[3][2])==floor(C[0][3]/C[3][3]) and floor(C[1][0]/C[3][0])==floor(C[1][1]/C[3][1])==floor(C[1][2]/C[3][2])==floor(C[1][3]/C[3][3]) and floor(C[2][0]/C[3][0])==floor(C[2][1]/C[3][1])==floor(C[2][2]/C[3][2])==floor(C[2][3]/C[3][3]) and sign(C[3][0])==sign(C[3][1])==sign(C[3][2])==sign(C[3][3]):
                qtrans1=floor(C[0][0]/C[3][0])
                qtrans2=floor(C[1][0]/C[3][0])
                qtrans3=floor(C[2][0]/C[3][0])
                exp[0].append(qtrans1)
                exp[1].append(qtrans2)
                exp[2].append(qtrans3)

                total_outputs+=1                           #total number of outputs
                inputs_array.append(input_counter)         #inputs required for this step
                total_inputs_array.append(total_inputs)    #total inputs required until this output
                input_counter=0                            

                trans_mat=matrix(RR,[[0,0,0,1],[1,0,0,-qtrans1],[0,1,0,-qtrans2],[0,0,1,-qtrans3]])
                C=trans_mat*C
                
            else:

                q1=randint(1,pq_lim)
                q2=randint(0,q1)
                q3=randint(0,q2)

                if i==0:
                    q2=randint(0,pq_lim)
                    q3=randint(0,pq_lim)
                
                starting_cf[0].append(q1)
                starting_cf[1].append(q2)
                starting_cf[2].append(q3)
  
                input_counter+=1
                total_inputs+=1

                trans_mat=matrix(RR,[[q1,1,0,0],[q2,0,1,0],[q3,0,0,1],[1,0,0,0]])
                C=C*trans_mat

        else:

                q1=randint(1,pq_lim)
                q2=randint(0,q1)
                q3=randint(0,q2)

                if i==0:
                    q2=randint(0,pq_lim)
                    q3=randint(0,pq_lim)
                
                starting_cf[0].append(q1)
                starting_cf[1].append(q2)
                starting_cf[2].append(q3)
  
                input_counter+=1
                total_inputs+=1

                trans_mat=matrix(RR,[[q1,1,0,0],[q2,0,1,0],[q3,0,0,1],[1,0,0,0]])
                C=C*trans_mat

        i+=1

    print("Total inputs:",total_inputs)
    print("Total inputs array:",total_inputs_array)
    print("Inputs array:",inputs_array)
    print("Total outputs:",total_outputs)

    print("Initial expansion is:",starting_cf)
    print("Final expansion:",exp)

    return total_inputs_array


def tri_plot_random_inputs_upto(out_lim,pq_lim,coeff_lim,N,col):

    inputs_ar=[]

    for i in range(N):
        c1=random_vector(4,coeff_lim)
        c2=random_vector(4,coeff_lim)
        c3=random_vector(4,coeff_lim)
        c4=random_vector(4,coeff_lim)
        C=matrix(RR,[c1,c2,c3,c4])

        while C.determinant()==0:
            c1=random_vector(4,coeff_lim)
            c2=random_vector(4,coeff_lim)
            c3=random_vector(4,coeff_lim)
            c4=random_vector(4,coeff_lim)
            C=matrix(RR,[c1,c2,c3,c4])

        total_inp_ar=tri_count_random_inputs_upto(C,out_lim,pq_lim)
        inputs_ar.append(total_inp_ar)

    mean_inputs_ar=[]

    for k in range(out_lim):
        mean=0
        for i in range(N):
            mean+=inputs_ar[i][k]
        
        mean=mean/N
        mean_inputs_ar.append(RR(mean))

    points_for_plot=[]
    for i in range(out_lim):
        points_for_plot.append([i+1,mean_inputs_ar[i]])
    print(out_lim)
    G=list_plot(points_for_plot, gridlines=True, plotjoined=True,color=col,thickness=1)
   
    return G


def tri_plot_random_inputs_without_mean_upto(out_lim,pq_lim,coeff_lim,N,col):

    inputs_ar=[]

    c1=random_vector(4,coeff_lim)
    c2=random_vector(4,coeff_lim)
    c3=random_vector(4,coeff_lim)
    c4=random_vector(4,coeff_lim)
    C=matrix(RR,[c1,c2,c3,c4])

    while C.determinant()==0:
        c1=random_vector(4,coeff_lim)
        c2=random_vector(4,coeff_lim)
        c3=random_vector(4,coeff_lim)
        c4=random_vector(4,coeff_lim)
        C=matrix(RR,[c1,c2,c3,c4])


    for i in range(N):
        total_inp_ar=tri_count_random_inputs_upto(C,out_lim,pq_lim)
        inputs_ar.append(total_inp_ar)

    
    seq_to_plot=[]

    for i in range(len(inputs_ar)):
        plot_to_app=[]
        for j in range(out_lim):
            plot_to_app.append([j+1,inputs_ar[i][j]])
        seq_to_plot.append(plot_to_app)

    G=list_plot(seq_to_plot[0], gridlines=True, plotjoined=True,color=col,thickness=1)

    for i in range(1,len(inputs_ar)):
        G=G+list_plot(seq_to_plot[i], gridlines=True, plotjoined=True,color=col,thickness=1)
    
    return G


#BILINEAR m=3

def tri_biquad_count_random_inputs_upto(C1,C2,C3,C4,out_lim,pq_lim):
    exp=[[],[],[]]

    starting_cf=[[],[],[]]

    [ix,iy]=[0,0]
    ninput=0

    x1=[]
    x2=[]
    x3=[]
    y1=[]
    y2=[]
    y3=[]

    for i in range(1000*out_lim):
        com1=randint(1,pq_lim)
        com2=randint(0,com1)
        com3=randint(0,com2)
        x1.append(com1)
        x2.append(com2)
        x3.append(com3)

    for i in range(1000*out_lim):
        com1=randint(1,pq_lim)
        com2=randint(0,com1)
        com3=randint(0,com2)
        y1.append(com1)
        y2.append(com2)
        y3.append(com3)

    i=0

    inputs_array=[]
    total_inputs_array=[]
    total_inputs=0
    total_outputs=0
    input_counter=0


    while total_outputs<out_lim and i<50:
        print("Step:",i)
        
        

        if C4[0][0]!=0 and C4[0][1]!=0 and C4[0][2]!=0 and C4[0][3]!=0 and C4[1][0]!=0 and C4[1][1]!=0 and C4[1][2]!=0 and C4[1][3]!=0 and C4[2][0]!=0 and C4[2][1]!=0 and C4[2][2]!=0 and C4[2][3]!=0 and C4[3][0]!=0 and C4[3][1]!=0 and C4[3][2]!=0 and C4[3][3]!=0:
            if floor(C1[0][0]/C4[0][0])==floor(C1[0][1]/C4[0][1])==floor(C1[0][2]/C4[0][2])==floor(C1[0][3]/C4[0][3])==floor(C1[1][0]/C4[1][0])==floor(C1[1][1]/C4[1][1])==floor(C1[1][2]/C4[1][2])==floor(C1[1][3]/C4[1][3])==floor(C1[2][0]/C4[2][0])==floor(C1[2][1]/C4[2][1])==floor(C1[2][2]/C4[2][2])==floor(C1[2][3]/C4[2][3])==floor(C1[3][0]/C4[3][0])==floor(C1[3][1]/C4[3][1])==floor(C1[3][2]/C4[3][2])==floor(C1[3][3]/C4[3][3]) and floor(C2[0][0]/C4[0][0])==floor(C2[0][1]/C4[0][1])==floor(C2[0][2]/C4[0][2])==floor(C2[0][3]/C4[0][3])==floor(C2[1][0]/C4[1][0])==floor(C2[1][1]/C4[1][1])==floor(C2[1][2]/C4[1][2])==floor(C2[1][3]/C4[1][3])==floor(C2[2][0]/C4[2][0])==floor(C2[2][1]/C4[2][1])==floor(C2[2][2]/C4[2][2])==floor(C2[2][3]/C4[2][3])==floor(C2[3][0]/C4[3][0])==floor(C2[3][1]/C4[3][1])==floor(C2[3][2]/C4[3][2])==floor(C2[3][3]/C4[3][3]) and floor(C3[0][0]/C4[0][0])==floor(C3[0][1]/C4[0][1])==floor(C3[0][2]/C4[0][2])==floor(C3[0][3]/C4[0][3])==floor(C3[1][0]/C4[1][0])==floor(C3[1][1]/C4[1][1])==floor(C3[1][2]/C4[1][2])==floor(C3[1][3]/C4[1][3])==floor(C3[2][0]/C4[2][0])==floor(C3[2][1]/C4[2][1])==floor(C3[2][2]/C4[2][2])==floor(C3[2][3]/C4[2][3])==floor(C3[3][0]/C4[3][0])==floor(C3[3][1]/C4[3][1])==floor(C3[3][2]/C4[3][2])==floor(C3[3][3]/C4[3][3]) and sign(C4[0][0])==sign(C4[0][1])==sign(C4[0][2])==sign(C4[0][3])==sign(C4[1][0])==sign(C4[1][1])==sign(C4[1][2])==sign(C4[1][3])==sign(C4[2][0])==sign(C4[2][1])==sign(C4[2][2])==sign(C4[2][3])==sign(C4[3][0])==sign(C4[3][1])==sign(C4[3][2])==sign(C4[3][3]):
                [qtrans1,qtrans2,qtrans3]=[floor(C1[0][0]/C4[0][0]),floor(C2[0][0]/C4[0][0]),floor(C3[0][0]/C4[0][0])]

                exp[0].append(qtrans1)
                exp[1].append(qtrans2)
                exp[2].append(qtrans3)

                total_outputs+=1                           #total number of outputs
                inputs_array.append(input_counter)         #inputs required for this step
                total_inputs_array.append(total_inputs)    #total inputs required until this output
                input_counter=0                            

               
                [C1,C2,C3,C4]=[C4,C1-qtrans1*C4,C2-qtrans2*C4,C3-qtrans3*C4]



            else:
                ninput+=1

                if ix==iy:
                    [q1,q2,q3]=[x1[ix],x2[ix],x3[ix]]
                    
                    trans_mat=matrix(RR,[[q1,q2,q3,1],[1,0,0,0],[0,1,0,0],[0,0,1,0]])


               
                    C1=trans_mat*C1
                    C2=trans_mat*C2
                    C3=trans_mat*C3
                    C4=trans_mat*C4

                    C1=C1.transpose()
                    C2=C2.transpose()
                    C3=C3.transpose()
                    C4=C4.transpose()

                    input_counter+=1
                    total_inputs+=1

                    ix+=1

                else:
                    [q1,q2,q3]=[y1[iy],y2[iy],y3[iy]]
                    
                    trans_mat=matrix(RR,[[q1,q2,q3,1],[1,0,0,0],[0,1,0,0],[0,0,1,0]])

                    C1=trans_mat*C1
                    C2=trans_mat*C2
                    C3=trans_mat*C3
                    C4=trans_mat*C4

                    C1=C1.transpose()
                    C2=C2.transpose()
                    C3=C3.transpose()
                    C4=C4.transpose()

                    input_counter+=1
                    total_inputs+=1

                    iy+=1
        else:
            ninput+=1

            if ix==iy:
                [q1,q2,q3]=[x1[ix],x2[ix],x3[ix]]
                
                trans_mat=matrix(RR,[[q1,q2,q3,1],[1,0,0,0],[0,1,0,0],[0,0,1,0]])

                C1=trans_mat*C1
                C2=trans_mat*C2
                C3=trans_mat*C3
                C4=trans_mat*C4

                C1=C1.transpose()
                C2=C2.transpose()
                C3=C3.transpose()
                C4=C4.transpose()


                input_counter+=1
                total_inputs+=1

                ix+=1

            else:
                [q1,q2,q3]=[y1[iy],y2[iy],y3[iy]]
                
                trans_mat=matrix(RR,[[q1,q2,q3,1],[1,0,0,0],[0,1,0,0],[0,0,1,0]])

                C1=trans_mat*C1
                C2=trans_mat*C2
                C3=trans_mat*C3
                C4=trans_mat*C4

                C1=C1.transpose()
                C2=C2.transpose()
                C3=C3.transpose()
                C4=C4.transpose()


                input_counter+=1
                total_inputs+=1

        i+=1

    print("Total inputs:",total_inputs)
    print("Total inputs array:",total_inputs_array)
    print("Inputs array:",inputs_array)
    print("Total outputs:",total_outputs)
    print("Final expansion:",exp)

    return total_inputs_array


def tri_biquad_plot_random_inputs_without_mean_upto(out_lim,pq_lim,coeff_lim,N,col):

    inputs_ar=[]

    c11=random_vector(4,coeff_lim)
    c12=random_vector(4,coeff_lim)
    c13=random_vector(4,coeff_lim)
    c14=random_vector(4,coeff_lim)

    c21=random_vector(4,coeff_lim)
    c22=random_vector(4,coeff_lim)
    c23=random_vector(4,coeff_lim)
    c24=random_vector(4,coeff_lim)

    c31=random_vector(4,coeff_lim)
    c32=random_vector(4,coeff_lim)
    c33=random_vector(4,coeff_lim)
    c34=random_vector(4,coeff_lim)

    c41=random_vector(4,coeff_lim)
    c42=random_vector(4,coeff_lim)
    c43=random_vector(4,coeff_lim)
    c44=random_vector(4,coeff_lim)

    C1=matrix(RR,[c11,c12,c13,c14])
    C2=matrix(RR,[c21,c22,c23,c24])
    C3=matrix(RR,[c31,c32,c33,c34])
    C4=matrix(RR,[c41,c42,c43,c44])
    

    while C1.determinant()==0:
        c11=random_vector(4,coeff_lim)
        c12=random_vector(4,coeff_lim)
        c13=random_vector(4,coeff_lim)
        c14=random_vector(4,coeff_lim)

        C1=matrix(RR,[c11,c12,c13,c14])

    while C2.determinant()==0:
        c21=random_vector(4,coeff_lim)
        c22=random_vector(4,coeff_lim)
        c23=random_vector(4,coeff_lim)
        c24=random_vector(4,coeff_lim)

        C2=matrix(RR,[c21,c22,c23,c24])

    while C3.determinant()==0:
        c31=random_vector(4,coeff_lim)
        c32=random_vector(4,coeff_lim)
        c33=random_vector(4,coeff_lim)
        c34=random_vector(4,coeff_lim)

        C3=matrix(RR,[c31,c32,c33,c34])

    while C4.determinant()==0:
        c41=random_vector(4,coeff_lim)
        c42=random_vector(4,coeff_lim)
        c43=random_vector(4,coeff_lim)
        c44=random_vector(4,coeff_lim)

        C4=matrix(RR,[c41,c42,c43,c44])


    for i in range(N):
        total_inp_ar=tri_biquad_count_random_inputs_upto(C1,C2,C3,C4,out_lim,pq_lim)
        inputs_ar.append(total_inp_ar)

    
    seq_to_plot=[]

    for i in range(len(inputs_ar)):
        plot_to_app=[]
        for j in range(out_lim):
            plot_to_app.append([j+1,inputs_ar[i][j]])
        seq_to_plot.append(plot_to_app)

    G=list_plot(seq_to_plot[0], gridlines=True, plotjoined=True,color=col,thickness=1)

    for i in range(1,len(inputs_ar)):
        G=G+list_plot(seq_to_plot[i], gridlines=True, plotjoined=True,color=col,thickness=1)
    
    return G

