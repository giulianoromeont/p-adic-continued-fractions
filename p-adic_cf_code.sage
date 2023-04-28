#PRIMI 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97

p=5
K = Qp(p, prec = 100, type = 'capped-rel', print_mode = 'series')

def Ruban_floor(x):	
		v=x.valuation()
		b=0
		for i in range (v,1):
			b=b+(x[i]*p^i)
			
		return K(b)
		

def Ruban_cf(x,lim):	
		b=Ruban_floor(x)
		exp=[Rational(b)]
		if x-b==0:
			flag=0
		else:
			flag=1	
		
		i=0
		while flag==1  and i<lim:
			i=i+1
			x=1/(x-b)
			b=Ruban_floor(x)
			exp.append(Rational(b))
			if x-b==0:
				flag=0
			else:
				flag=1
	
		return exp
		
def Schneider_cf(x,lim):
	part_num=[]
	part_den=[]
	
	b=x[0]
	part_den.append(b)
	
	e=(x-b).valuation()
	a=p^e
	part_num.append(a)
	
	if x-b==0:
		flag=0
	else:
		flag=1	
	i=0
	while flag==1 and i<lim:
		i=i+1
		x=(p^e)/(x-b)
		b=x[0]
		part_den.append(b)
		e=(x-b).valuation()
		a=p^e
		part_num.append(a)
	return part_num,part_den
		
	


def Browkin_s(x):
		v=x.valuation()
		b=0
		t=x[v]
		if x[v]>(p-1)/2:
				t=-(p-x[v])
				carry=1
		else:
				carry=0
		b=b+(t*p^v)
		
		for i in range (v+1,1):
			t=x[i]+carry
			if t>(p-1)/2:
				t=-(p-t)
				carry=1
			else:
				carry=0
			b=b+(t*p^i)	
		return b
		
		
def Browkin_t(x):
		v=x.valuation()
		b=0
		t=x[v]
		if x[v]>(p-1)/2:
				t=-(p-x[v])
				carry=1
		else:
				carry=0
		b=b+(t*p^v)
		
		for i in range (v+1,0):
			t=x[i]+carry
			if t>(p-1)/2:
				t=-(p-t)
				carry=1
			else:
				carry=0
			b=b+(t*p^i)	
		return b
	
	
def BrowkinI_cf(x,lim):	
		b=Browkin_s(x)
		exp=[Rational(b)]
		if x-b==0:
			flag=0
		else:
			flag=1	
		
		i=0
		while flag==1  and i<lim:
			i=i+1
			x=1/(x-b)
			b=Browkin_s(x)
			exp.append(Rational(b))
			if x-b==0:
				flag=0
			else:
				flag=1
	
		return exp
	
	
def check_sign(x):	
		v=x.valuation()
		t=x[v]
		if x[v]>(p-1)/2:
				t=-(p-x[v])
				carry=1
		else:
				carry=0
		
		for i in range (v+1,1):
			t=x[i]+carry
			if t>(p-1)/2:
				t=-(p-t)
				carry=1
			else:
				carry=0
		if t==0:
			a0=1
		else:
			a0=0
		return a0

		
def BrowkinII_cf_Rat(x,lim):	
		b=Browkin_s(x)
		exp=[Rational(b)]
		if x-b==0:
			flag=0
		else:
			flag=1	
		
		i=0
		while flag==1 and i<lim:
			x=1/(x-b)
			i=i+1
			if i%2==0:
				b=Browkin_s(x)
			else:
				q=Browkin_t(x)
				b=q-check_sign(K(x))*sign(q)	
			exp.append(Rational(b))		
			if x-b==0:
				flag=0
			else:
				flag=1
						
		return exp

		
def New_cf_Rat(x,lim):	
		b=Browkin_s(x)
		exp=[Rational(b)]
		if x-b==0:
			flag=0
		else:
			flag=1	
		
		i=0
		while flag==1 and i<lim:
			i=i+1
			x=1/(x-b)
			if i%2==0:
				b=Browkin_s(x)
			else:
				b=Browkin_t(x)
					
			exp.append(Rational(b))		
			if x-b==0:
				flag=0
			else:
				flag=1					
		return exp


def BrowkinI_cf_Sqrt(x,lim):
		if K(x).is_square()==1:
			x=K(x).square_root()
			b=Browkin_s(x)
			exp=[Rational(b)]
			if x-b==0:
				flag=0
			else:
				flag=1	
			
			i=0
			while flag==1 and i<lim:
				x=1/(x-b)
				i=i+1
				b=Browkin_s(x)	
				exp.append(Rational(b))		
				if x-b==0:
					flag=0
				else:
					flag=1
		else:
			return print(x, 'is not a quadratic residue modulo p')
			
		return exp


def BrowkinII_cf_Sqrt(x,lim):
		if K(x).is_square()==1:
			x=K(x).square_root()
			b=Browkin_s(x)
			exp=[Rational(b)]
			if x-b==0:
				flag=0
			else:
				flag=1	
			
			i=0
			while flag==1 and i<lim:
				x=1/(x-b)
				i=i+1
				if i%2==0:
					b=Browkin_s(x)
				else:
					q=Browkin_t(x)
					b=q-check_sign(K(x))*sign(q)	
				exp.append(Rational(b))		
				if x-b==0:
					flag=0
				else:
					flag=1
		else:
			return print(x, 'is not a quadratic residue modulo p')
			
		return exp
		

def New_cf_Sqrt(x,lim):
		if K(x).is_square()==1:
			x=K(x).square_root()
			b=Browkin_s(x)
			exp=[Rational(b)]
			if x-b==0:
				flag=0
			else:
				flag=1	
			
			i=0
			while flag==1 and i<lim:
				x=1/(x-b)
				i=i+1
				if i%2==0:
					b=Browkin_s(x)
				else:
					b=Browkin_t(x)
				
				exp.append(Rational(b))		
				if x-b==0:
					flag=0
				else:
					flag=1
		else:
			return print(x, 'is not a quadratic residue modulo p')
			
		return exp
		
		

def detect_period_BrowkinI_Sqrt(D,lim):
		if K(D).is_square()==1:
			x=K(D).square_root()
			x1=x
			b1=Browkin_s(x)
			exp1=[Rational(b1)]
			
			P1=b1
			Q1=D-b1^2
			
			L1=[[0,1],[P1,Q1]]
			
			i=0
			
			flag1=0
			
			while flag1==0 and i<lim:
				x1=1/(x1-b1)
				
				i=i+1
				b1=Browkin_s(x1)
	
				exp1.append(Rational(b1))
				
				P1=b1*Q1-P1
				Q1=(D-P1^2)/Q1
				L1.append([P1,Q1])
				
				l=len(L1)
				j=0
				
				while flag1==0 and j<l-1:
					if L1[j]==[P1,Q1]:
						flag1=1
						length_period1=l-j-1
						length_preperiod1=j
					
					j=j+1
					
			if flag1==0:
				length_preperiod1=0
				length_period1=0		
				
			return [exp1,length_preperiod1,length_period1]
				
		else:
			return print(x, 'is not a quadratic residue modulo p')
		
		

def detect_period_BrowkinII_Sqrt(D,lim):
		if K(D).is_square()==1:
			x=K(D).square_root()
			x2=x
			b2=Browkin_s(x)
			exp2=[Rational(b2)]
			
			P2=b2
			Q2=D-b2^2
			L2=[[0,1],[P2,Q2]]
			
			i=0
			flag2=0
			
			while flag2==0 and i<lim:
				x2=1/(x2-b2)
				
				i=i+1
				
				if i%2==0:
					b2=Browkin_s(x2)
				else:
					q2=Browkin_t(x2)
					b2=q2-check_sign(K(x2))*sign(q2)
					
				exp2.append(Rational(b2))
				
				P2=b2*Q2-P2
				Q2=(D-P2^2)/Q2
				L2.append([P2,Q2])
				
				l=len(L2)
				j=0
				
				while flag2==0 and j<l-1:
					if L2[j]==[P2,Q2]:
						flag2=1
						length_period2=l-j-1
						length_preperiod2=j
					
					j=j+1
			
			if flag2==0:
				length_preperiod2=0
				length_period2=0			
				
			return [exp2,length_preperiod2,length_period2]
				
		else:
			return print(x, 'is not a quadratic residue modulo p')




def detect_period_New_Sqrt(D,lim):
		if K(D).is_square()==1:
			x=K(D).square_root()
			x3=x
			b3=Browkin_s(x3)
			exp3=[Rational(b3)]
			
			P3=b3
			Q3=D-b3^2
			L3=[[0,1],[P3,Q3]]
			
			i=0
			flag3=0
			
			while flag3==0 and i<lim:
				x3=1/(x3-b3)
				
				i=i+1
				
				if i%2==0:
					b3=Browkin_s(x3)
				else:
					b3=Browkin_t(x3)
					
				exp3.append(Rational(b3))
				
				P3=b3*Q3-P3
				Q3=(D-P3^2)/Q3
				L3.append([P3,Q3])
				
				l=len(L3)
				j=0
				
				while flag3==0 and j<l-1:
					if L3[j]==[P3,Q3]:
						flag3=1
						length_period3=l-j-1
						length_preperiod3=j
					
					j=j+1
			
			if flag3==0:
				length_preperiod3=0
				length_period3=0		
				
			return [exp3,length_preperiod3,length_period3]
				
		else:
			return print(x, 'is not a quadratic residue modulo p')



def try_roots_BrI(roof,lim):
	roots=[]
	for i in range(1,roof+1):
		if K(i).is_square()==True and sqrt(i)!=floor(sqrt(i)) and K(i).valuation()==0:
			roots.append(i)
	
	length_preperiods=[]
	length_periods=[]
	periodic_roots=[]
	number_periodic=0
	
	for i in range(len(roots)):
		x=detect_period_BrowkinI_Sqrt(roots[i],lim)
		if [x[1],x[2]]!=[0,0]:
			length_preperiods.append(x[1])
			length_periods.append(x[2])
			periodic_roots.append(roots[i])
			number_periodic=number_periodic+1
	
	return [length_preperiods,length_periods,periodic_roots,number_periodic]
	

def try_roots_BrII(roof,lim):
	roots=[]
	for i in range(1,roof+1):
		if K(i).is_square()==True and sqrt(i)!=floor(sqrt(i)) and K(i).valuation()==0:
			roots.append(i)
	
	length_preperiods=[]
	length_periods=[]
	periodic_roots=[]
	number_periodic=0
	
	for i in range(len(roots)):
		x=detect_period_BrowkinII_Sqrt(roots[i],lim)
		if [x[1],x[2]]!=[0,0]:
			length_preperiods.append(x[1])
			length_periods.append(x[2])
			periodic_roots.append(roots[i])
			number_periodic=number_periodic+1
	
	return [length_preperiods,length_periods,periodic_roots,number_periodic]


def try_roots_New(roof,lim):
	roots=[]
	for i in range(1,roof+1):
		if K(i).is_square()==True and sqrt(i)!=floor(sqrt(i)) and K(i).valuation()==0:
			roots.append(i)
	
	length_preperiods=[]
	length_periods=[]
	periodic_roots=[]
	number_periodic=0
	
	for i in range(len(roots)):
		x=detect_period_New_Sqrt(roots[i],lim)
		if [x[1],x[2]]!=[0,0] and x[2]%2==0:
			length_preperiods.append(x[1])
			length_periods.append(x[2])
			periodic_roots.append(roots[i])
			number_periodic=number_periodic+1
	
	return [length_preperiods,length_periods,periodic_roots,number_periodic]
	
def count_periods(roof,lim):

	total=0
	for i in range(1,roof+1):
		if K(i).is_square()==True and sqrt(i)!=floor(sqrt(i)) and K(i).valuation()==0:
			total=total+1
	
	X=try_roots_BrI(roof,lim)
	Y=try_roots_BrII(roof,lim)
	Z=try_roots_New(roof,lim)
	
	mean_preperiod_Br1=mean(X[0])
	mean_preperiod_Br2=mean(Y[0])
	mean_preperiod_New=mean(Z[0])
	
	mean_period_Br1=mean(X[1])
	mean_period_Br2=mean(Y[1])
	mean_period_New=mean(Z[1])
	
	count_periodic_Br1=X[3]
	count_periodic_Br2=Y[3]
	count_periodic_New=Z[3]
	
	
	return 'PRIME:', p, 'PREPERIODS:',[float(mean_preperiod_Br1),float(mean_preperiod_Br2),float(mean_preperiod_New)],'PERIODS:',[float(mean_period_Br1),float(mean_period_Br2),float(mean_period_New)],'PERIODIC:',[count_periodic_Br1,count_periodic_Br2,count_periodic_New],'PERIODS LIST:', [X[1],Y[1],Z[1]],'PERIODIC ROOTS LIST:',[X[2],Y[2],Z[2]] ,'TOTAL:',total
	

def approximation_Sqrt(D,lim):
	A1=BrowkinI_cf_Sqrt(D,lim)
	A2=BrowkinII_cf_Sqrt(D,lim)
	A3=New_cf_Sqrt(D,lim)
	
	v1=0
	v2=0
	v3=0
	
	for i in range(lim):
		r1=A1[i]
		r2=A2[i]
		r3=A3[i]
		v1=v1+K(r1).valuation()
		v2=v2+K(r2).valuation()
		v3=v3+K(r3).valuation()
		
	return [v1,v2,v3]
	
	
def quality_approximation(roof):
	roots=[]
	for i in range(1,roof+1):
		if K(i).is_square()==True and sqrt(i)!=floor(sqrt(i)) and K(i).valuation()==0:
			roots.append(i)

	BI=[[],[],[]]
	BII=[[],[],[]]
	New=[[],[],[]]

	for i in range(len(roots)):
		app10=approximation_Sqrt(roots[i],10)
		BI[0].append(float(app10[0]))
		BII[0].append(float(app10[1]))
		New[0].append(float(app10[2]))

		app100=approximation_Sqrt(roots[i],100)
		BI[1].append(float(app100[0]))
		BII[1].append(float(app100[1]))
		New[1].append(float(app100[2]))

		app1000=approximation_Sqrt(roots[i],1000)
		BI[2].append(float(app1000[0]))
		BII[2].append(float(app1000[1]))
		New[2].append(float(app1000[2]))
		
	mean1=[mean(BI[0]),mean(BI[1]),mean(BI[2])]
	mean2=[mean(BII[0]),mean(BII[1]),mean(BII[2])]
	mean3=[mean(New[0]),mean(New[1]),mean(New[2])]
	
	return mean1,mean2,mean3