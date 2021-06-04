# f1,..,f28 forms a basis of I_2: Since the degree (in t_i) is 2, this is enough. Symmetric bilinear forms are defined on e_i\times e_j.
# 
# T1 is the matrix containing the first rows of all matrices M1,...,M6.
# T2 contains the first columns.
# T3 contains the second rows.
# T4 contains the second columns.
# 
# Mini contains "all" bordered determinants.
# 
# We generate J by all these bordered determinants, I_2 and the products of 2x2 minors of the T_i.
# These are not needed to generate the radical but are needed for the Groebner basis.
# They are contained in the radical by looking at the Fano scheme of P1 x P1.
                                             
sage: R.<x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4,w1,w2,w3,w4,v1,v2,v3,v4,c1,c2,c3,c4,b1,b2,b3,b4,t1,t2,t3,t4,t5,t6,t7>=PolynomialRing(QQ)
sage: M1=matrix([[x1,x3],[x2,x4]])
sage: M2=matrix([[y1,y3],[y2,y4]])
sage: M3=matrix([[z1,z3],[z2,z4]])
sage: M4=matrix([[w1,w3],[w2,w4]])
sage: M5=matrix([[v1,v3],[v2,v4]])
sage: M6=matrix([[c1,c3],[c2,c4]])
sage: M7=matrix([[b1,b3],[b2,b4]])
sage: M=t1*M1+t2*M2+t3*M3+t4*M4+t5*M5+t6*M6+t7*M7
sage: d=det(M)
sage: f1=d.subs(t1=1,t2=0,t3=0,t4=0,t5=0,t6=0,t7=0)
sage: f7=d.subs(t1=0,t2=1,t3=0,t4=0,t5=0,t6=0,t7=0)
sage: f12=d.subs(t1=0,t2=0,t3=1,t4=0,t5=0,t6=0,t7=0)
sage: f16=d.subs(t1=0,t2=0,t3=0,t4=1,t5=0,t6=0,t7=0)
sage: f19=d.subs(t1=0,t2=0,t3=0,t4=0,t5=1,t6=0,t7=0)
sage: f21=d.subs(t1=0,t2=0,t3=0,t4=0,t5=0,t6=1,t7=0)
sage: f2=d.subs(t1=1,t2=1,t3=0,t4=0,t5=0,t6=0,t7=0)-f1-f7
sage: f3=d.subs(t1=1,t2=0,t3=1,t4=0,t5=0,t6=0,t7=0)-f1-f12
sage: f4=d.subs(t1=1,t2=0,t3=0,t4=1,t5=0,t6=0,t7=0)-f1-f16
sage: f5=d.subs(t1=1,t2=0,t3=0,t4=0,t5=1,t6=0,t7=0)-f1-f19
sage: f6=d.subs(t1=1,t2=0,t3=0,t4=0,t5=0,t6=1,t7=0)-f1-f21
sage: f8=d.subs(t1=0,t2=1,t3=1,t4=0,t5=0,t6=0,t7=0)-f7-f12
sage: f9=d.subs(t1=0,t2=1,t3=0,t4=1,t5=0,t6=0,t7=0)-f7-f16
sage: f10=d.subs(t1=0,t2=1,t3=0,t4=0,t5=1,t6=0,t7=0)-f7-f19
sage: f11=d.subs(t1=0,t2=1,t3=0,t4=0,t5=0,t6=1,t7=0)-f7-f21
sage: f13=d.subs(t1=0,t2=0,t3=1,t4=1,t5=0,t6=0,t7=0)-f12-f16
sage: f14=d.subs(t1=0,t2=0,t3=1,t4=0,t5=1,t6=0,t7=0)-f12-f19
sage: f15=d.subs(t1=0,t2=0,t3=1,t4=0,t5=0,t6=1,t7=0)-f12-f21
sage: f17=d.subs(t1=0,t2=0,t3=0,t4=1,t5=1,t6=0,t7=0)-f16-f19
sage: f18=d.subs(t1=0,t2=0,t3=0,t4=1,t5=0,t6=1,t7=0)-f16-f21
sage: f20=d.subs(t1=0,t2=0,t3=0,t4=0,t5=1,t6=1,t7=0)-f19-f21
sage: f22=d.subs(t1=0,t2=0,t3=0,t4=0,t5=0,t6=0,t7=1)
sage: f23=d.subs(t1=0,t2=0,t3=0,t4=0,t5=0,t6=1,t7=1)-f22-f21
sage: f24=d.subs(t1=0,t2=0,t3=0,t4=0,t5=1,t6=0,t7=1)-f22-f19
sage: f25=d.subs(t1=0,t2=0,t3=0,t4=1,t5=0,t6=0,t7=1)-f22-f16
sage: f26=d.subs(t1=0,t2=0,t3=1,t4=0,t5=0,t6=0,t7=1)-f22-f12
sage: f27=d.subs(t1=0,t2=1,t3=0,t4=0,t5=0,t6=0,t7=1)-f22-f7
sage: f28=d.subs(t1=1,t2=0,t3=0,t4=0,t5=0,t6=0,t7=1)-f22-f1
sage: I=ideal(f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,f19,f20,f21,f22,f23,f24,f25,f26,f27,f28)
sage: T=matrix([[x1,x2,x3,x4],[y1,y2,y3,y4],[z1,z2,z3,z4],[w1,w2,w3,w4],[v1,v2,v3,v4],[c1,c2,c3,c4],[b1,b2,b3,b4]])
sage: T1=T[[0..6],[0,1]]
sage: T2=T[[0..6],[0,2]]
sage: T3=T[[0..6],[2,3]]
sage: T4=T[[0..6],[1,3]]
sage: IT1=ideal(T1.minors(2))
sage: IT2=ideal(T2.minors(2))
sage: IT3=ideal(T3.minors(2))
sage: IT4=ideal(T4.minors(2))

sage: Z=[]
sage: Z.append(M1)
sage: Z.append(M2)
sage: Z.append(M3)
sage: Z.append(M4)
sage: Z.append(M5)
sage: Z.append(M6)
sage: Z.append(M7)

sage: mini=[]
sage: for i in [0..6]:   
....:     for j in [0..6]:    
....:         for k in [0..6]:
....:             Q=x1*zero_matrix(3)
....:             Q[[0..1],[0..1]]=Z[i]
....:             Q[[0..1],2]=Z[j][[0..1],0]
....:             Q[2,[0..1]]=Z[k][0,[0..1]]
....:             mini.append(Q.det())
....:             Q[[0..1],2]=Z[j][[0..1],0]
....:             Q[2,[0..1]]=Z[k][1,[0..1]]
....:             mini.append(Q.det())
....:             Q[[0..1],2]=Z[j][[0..1],1]
....:             Q[2,[0..1]]=Z[k][0,[0..1]]
....:             mini.append(Q.det())
....:             Q[[0..1],2]=Z[j][[0..1],1]
....:             Q[2,[0..1]]=Z[k][1,[0..1]]
....:             mini.append(Q.det())
sage: J=I+R.ideal(mini)+(IT1*IT4)+(IT2*IT3)
sage: J.basis_is_groebner()
