class ComputingAs_TwoBridgeLink:
    def InverseMod(self, q, p):  # compute the inverse of q modulus p, provided p,q are coprime
        k=0
        while (k*q)%p!=1:
            k+=1
        return k 

    def ListAdd(self, l1, l2):  # vector addition
        for i in range(len(l1)):
            l1[i]+=l2[i]
        return l1

    def SymAlexanderGrading(self, p, q):    # compute the symmetric Alexander grading of link Floer complex of the two-bridge link
        A=[[0,0] for i in range(2*p)]
        l=self.InverseMod(q,2*p)
        k=[0,0]
        for i in range(2*p):
            if i%2==0:
                k=self.ListAdd(k,[(-1)**((-i*q)//p),0]) 
            else:
                k=self.ListAdd(k,[0,(-1)**(-(i*q+1)//p)])
            A[i][0]=k[0]
            A[i][1]=k[1]
        R1=[A[i][0] for i in range(2*p)]
        R2=[A[i][1] for i in range(2*p)]
        k=[(min(R1)+max(R1))/2, (min(R2)+max(R2))/2]
        for i in range(2*p):
            for j in range(2):
                A[i][j]-=k[j]
        return A
    
    def NetMesh(self, X):     # compute the mesh of a list
        S=sorted(X)
        U=[S[i+1]-S[i] for i in range(len(S)-1)]
        return min(U)

    def NetMesh2(self, X, Y, p):
        if X==[]:
            if len(Y)==1:
                return p-max(Y)-1
            else:
                return min(self.NetMesh(Y)//2, p-max(Y)-1)
        elif Y==[]:
            if len(X)==1:
                return min(X[0]-1, p-X[0]-1)
            else:
                return min(self.NetMesh(X)//2, min(X)-1, p-max(X)-1)
        else:
            if len(X)==1 and len(Y)==1:
                return min(X[0]-1, p-X[0]-1, p-Y[0]-1)
            elif len(X)==1 and len(Y)>1:
                return min([min(X)-1, p-max(X)-1, self.NetMesh(Y)//2, p-max(Y)-1])
            elif len(X)>1 and len(Y)==1:
                return min([self.NetMesh(X)//2, min(X)-1, p-max(X)-1, p-max(Y)-1])
            else:
                return min([self.NetMesh(X)//2, min(X)-1, p-max(X)-1, self.NetMesh(Y)//2, p-max(Y)-1])

    def BX(self, i, j, p, q):
        return [(i+j)%(2*p), (i-1-j)%(2*p)]

    def BY(self, i, j, p, q):
        return [(q-i-1-j)%(2*p), (q-i+j)%(2*p)]
    
    def Mod(self, m, n, d):
        # define the remainder of division m by n by offset d
        # d might be negative
        res=m%n
        while res>=d+n:
            res-=n
        while res<d:
            res+=n
        return res
        

    def Differential(self, p, q): # getting the differential of link Floer chain complex
        D={}
        for i in range(2*p):
            D[i]=set()
        
        X=[]
        Y=[0]
        for i in range(p):
            k=self.NetMesh2(X,Y,p)
            for j in range(k+1):
                Bigon=[((-1)**i*(i + ((-1)**i + 1)/2)*q + (-1)**(i + 1)*j - ((-1)**i + 1)/2)%(2*p), ((-1)**i*(i + ((-1)**i + 1)/2)*q + (-1)**i*j + ((-1)**i - 1)/
 2)%(2*p)]
                if Bigon[1] in D[Bigon[0]]:
                    D[Bigon[0]].remove(Bigon[1])
                else:
                    D[Bigon[0]].add(Bigon[1])
            if i%2==0:
                X.append(abs(self.Mod((i + 1)*q, 2*p, -p + 1)))
            else:
                Y.append(abs(self.Mod((i + 1)*q, 2*p, -p + 1)))
        X=[0]
        Y=[]
        for i in range(p):
            #print 'i=', i, 'X=', X, ',Y=', Y
            k=self.NetMesh2(Y,X,p)
            if i%2==0:
                for j in range(k+1):
                    Bigon=self.BX(-i*q, j, p, q)
                    if Bigon[1] in D[Bigon[0]]:
                        D[Bigon[0]].remove(Bigon[1])
                    else:
                        D[Bigon[0]].add(Bigon[1])
            else:
                for j in range(k+1):
                    Bigon=self.BY(-i*q, j, p, q)
                    if Bigon[1] in D[Bigon[0]]:
                        D[Bigon[0]].remove(Bigon[1])
                    else:
                        D[Bigon[0]].add(Bigon[1])                    
            if i%2==1:
                X.append(abs(self.Mod((i + 1)*q, 2*p, -p + 1)))
            else:
                Y.append(abs(self.Mod((i + 1)*q, 2*p, -p + 1)))
        X=[]
        for key in D:
            if D[key]!=[]:
                X.append(key)
        for i in X:
            for j in D[i]:
                D[(i+p)%(2*p)].add((j+p)%(2*p))
        return D

    def AsHat(self, p, q, s1, s2):
        # compute the generalized Link Floer complex for hat version
        SA=self.SymAlexanderGrading(p,q)
        grading={}
        for i in range(2*p):
            grading[i]=SA[i][0]+SA[i][1]-2*max(SA[i][0]-s1,0)-2*max(SA[i][1]-s2,0)
        D=self.Differential(p,q)
        for i in D:
            X=set()
            for j in D[i]:
                k=(grading[j]-grading[i]+1)/2
                if k>0:
                    X.add(j)
            D[i]=D[i]-X
        return D

    
                
            
