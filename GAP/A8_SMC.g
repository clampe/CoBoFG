LoadPackage("QPA");

p:=3;
k:=FiniteField(p*p);

#Set up the path algebra
Q:=Quiver(5, [ [1,5,"a1"], [2,5,"a2"], [3,5,"a3"], [4,5,"a4"], [5,1,"b1"], [5,2,"b2"], [5,3,"b3"], [5,4,"b4"] ]);
kQ:=PathAlgebra(k, Q);
AssignGeneratorVariables(kQ);
relations:=[ a1*b3, a2*b4, a3*b1, a4*b2, a1*b1*a1-a1*b4*a4, a2*b3*a3-a2*b2*a2 , a3*b3*a3-a3*b4*a4 , a4*b3*a3-a4*b4*a4, b1*a1*b1-b4*a4*b1, b3*a3*b2-b2*a2*b2, b3*a3*b3-b4*a4*b3, b3*a3*b4-b4*a4*b4, b4*a4+b2*a2+b1*a1+b3*a3-b1*a1*b1*a1 ];
gb:=GBNPGroebnerBasis(relations,kQ);
I:=Ideal(kQ, gb);
GroebnerBasis(I,gb);
B:=kQ/I;
Print("Dimension ",Dimension(B),"\n");

#Calculate simple and projective modules
irr:=SimpleModules(B);
l:=Length(irr);
#We assume: Thesis <-> Calculation
# 		 0 <-> [1,0,0,0,0]
# 		 1 <-> [0,1,0,0,0]
# 		 2 <-> [0,0,1,0,0]
# 		 3 <-> [0,0,0,1,0]
# 		 4 <-> [0,0,0,0,1]
for M in irr do
	Print("Dimension of simple:", Dimension(M), "\n");
od;
proj:=IndecProjectiveModules(B);
Print("\nRadical Series of projectives:\n");
for M in proj do
	Print(RadicalSeries(M)," dim",Dimension(M),"\n");
od;

#Find the restrictions
res:=[0,0,0,0,0];
coordinates:=[ [1,1,6], [3,3,6], [2,2,4], [4,4,6], [5,5,12] ];
n:=[0,0,1,0,0];
#Morphisms between Projectives
for i in [1..l] do
	for j in [1..l] do
		#Print("\nHom(",i,",",j,")\n");
		for morphism in HomOverAlgebra(proj[i], proj[j]) do
			M:=Kernel(morphism);
			pos:=Positions(coordinates, [i,j,Dimension(M)]);
			if pos<>[] then res[pos[1]]:=Image(morphism);fi;
			#Print(RadicalSeries(Image(morphism)),"\n");
		od;
	od;
od;

Print("\nRestrictions\n");
for i in [1..l] do
	Print("S_",i-1,"\n");
	Print(RadicalSeries(res[i]),"\n");
	Print(SocleSeries(res[i]),"\n");
od;

#Compile the simple minded collection
SMC:=[];
Print("\nCalculating Szyzygies\n");
for i in [1..l] do
	M:=res[i];	
	for j in [1..n[i]] do
		M:=1stSyzygy(M);
	od;
	Add(SMC, M);
od;
Print("\nSimple minded collection by their radical and socle series\n");
for i in [1..l] do
	M:=SMC[i];
	Print("X_",i-1,"\n");
	Print(RadicalSeries(M),"\n");
	Print(SocleSeries(M),"\n");
od;

#Check the conditions posed by the theorem on the SMC
Print("\nChecking conditions (i) and (ii)\n");
pass:=true;
Print("(i) ");
for i in [1..l] do
	M:=SMC[i];
	if Length(HomOverAlgebra(M,M))<> 1 then pass:=fail;fi;
	#Print(Length(HomOverAlgebra(M,M)),"\n");
od;
if pass=true then Print("true\n");fi;
pass:=true;
Print("(ii) ");
for i in [1..l] do
	for j in [1..l] do
		M:=SMC[i];
		N:=SMC[j];
		if i<>j and n[i]<=n[j] then
			if Length(HomOverAlgebra(M,N))<> 0 then pass:=fail;fi;
			#Print(i," ",j," dim ",Length(HomOverAlgebra(M,N)),"\n");
		fi;
	od;
od;
if pass=true then Print("true\n");fi;
QUIT_GAP();
