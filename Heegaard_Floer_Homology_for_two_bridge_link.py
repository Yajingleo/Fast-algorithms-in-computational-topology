class ComputingAs_TwoBridgeLink:
    def InverseMod(q, p):  # compute the inverse of q modulus p, provided p,q are coprime
	k=0
	while (k*q)%p!=1:
		k+=1
	return k
    
    def Delta(i, p, q):    # compute the 
	if i%2==0:
		return 
	else:
		return 

    def ListAdd(l1, l2):    # vector addition
	for i in range(len(l1)):
		l1[i]+=l2[i]
	return l1

    def SATable(p, q):      # compute the symmetric Alexander grading of link Floer complex of the two-bridge link 
	A=[[0,0] for i in range(2*p)]
	l=InverseMod(q,2*p)
	k=[0,0]
	for i in range(2*p):
            if i%2==0:
		k=ListAdd(k,[(-1)**((-i*q)//p),0])
	    else:
                k=ListAdd(k,[0,(-1)**(-(i*q+1)//p)])
            A[i][0]=k[0]
	    A[i][1]=k[1]
	R1=[A[i][0] for i in range(2*p)]
	R2=[A[i][1] for i in range(2*p)]
	k=[(min(R1)+max(R1))/2, (min(R2)+max(R2))/2]
	for i in range(2*p):
		for j in range(2):
			A[i][j]-=k[j]
	return A
    
    def NetMesh(X):     # compute the mesh of a list
	S=sorted(X)
	U=[S[i+1]-S[i] for i in range(len(S)-1)]
	return min(U)

    def NetMesh2(X, Y, p): 
	return min([NetMesh(X)//2, min(X)-1, p-max(X)-1, NetMesh(Y)//2, p-max(Y)-1])

    def BX(i, j, p, q):
	return [(i+j)%(2*p), (i-1-j)%(2*p)]

    def BY(i, j, p, q):
	return [(q-i-1-j)%(2*p), (q-i+j)%(2*p)]
    
    def Mod(m, n, d):   # define the remainder of division m by n by offset d
        if (m%n>=d):
            return m%n
        else:
            return m%n+n
        

    def Differetial(p, q): # getting the differential of link Floer chain complex
        D={}
        X=[]
        Y=[0]
        for i in range(p):
            k=NetMesh(X,Y,p)
            for j in range(k+1):
                Bigon=[((-1)**i*(i + ((-1)**i + 1)/2)*q + (-1)**(i + 1)*j - ((-1)**i + 1)/2)//(2*p), ((-1)**i*(i + ((-1)**i + 1)/2)*q + (-1)**i*j + ((-1)**i - 1)/
 2)//(2*p)]
                if Bigon[1] in D[Bigon[0]]:
                    D[Bigon[0]].remove(Bigon[1])
                else:
                    D[Bigon[0]].add(Bigon[1])
                if i%2==0:
                    X.append(abs(Mod((i + 1)*q, 2*p, -p + 1)))
                else:
                    Y.append(abs(Mod((i + 1)*q, 2*p, -p + 1)))
        X=[0]
        Y=[]
        for i in range(p):
            k=NetMesh2(Y,X,p)
            if i%2==0:
                for j in range(k+1):
                    Bigon=BX(-i*q, j, p, q)
                    if Bigon[1] in D[Bigon[0]]:
                        D[Bigon[0]].remove(Bigon[1])
                    else:
                        D[Bigon[0]].add(Bigon[1])
            else:
                for j in range(k+1):
                    Bigon=BY(-i*q, j, p, q)
                    if i%2==1:
                        X.append(abs(Mod((i + 1)*q, 2*p, -p + 1)))
                    else:
                        Y.append(abs(Mod((i + 1)*q, 2*p, -p + 1)))
        for key in D:
            for i in D[key]:
                D[(key+p)//(2*p)].add((i+p)//(2*p))
        return D

        
        

        
                
                
                    
        
        
