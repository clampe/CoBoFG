#Von BS
CenterOfGroupAlgebra:=function(G,p) #Z(F_pG) where F_p is the field with q elements
local ct,dim,i,j,l,k,SCT;
	if IsCharacterTable(G) then ct:=G; else ct:=CharacterTable(G); Irr(ct); fi;
	dim:=NrConjugacyClasses(ct);
	SCT:=EmptySCTable(dim, Zero(GF(p)), "symmetric");
	for i in [1..dim] do
		for j in [i..dim] do
			l:=[];
			for k in [1..2*dim] do
				if k mod 2=1 then l[k]:=ClassMultiplicationCoefficient(ct,i,j,(k+1)/2)*One(GF(p)); else l[k]:=k/2; fi;
			od;
			SetEntrySCTable(SCT,i,j,l);
		od;
	od;
	return AlgebraByStructureConstants(GF(p),SCT);
end;

#von BS
#Calculates the radical/Loewy layers of an algebra
LoewyLayers:=function(A)
	local J,JJ,layers,last;
	J:=RadicalOfAlgebra(A);
	JJ:=J;
	last:=A;
	layers:=[A/J];
	while Dimension(JJ)>0 do
		last:=JJ;
		JJ:=ProductSpace(J,JJ);
		Add(layers, last/JJ);
	od;
	return layers;
end;

DimLoewyLayers:=function(A)
	local layer, dim;
	dim:=[];
	for layer in LoewyLayers(A) do
		Add(dim, Dimension(layer));
	od;
	return dim;
end;

#Lists all normal p-blocks using the method described in the notes. Also calculates k(B), l(B) and the dimension of LoewyLayers of the center
#A first version of this is due to B. Sambale.
ListNormalBlocks:=function(p)
	local SG, F, H, G, S, ct, pb, mct, i, sc, Q, Z, dim, ci, cis;
	SG:=List(ConjugacyClassesSubgroups(GL(2,p)),Representative);;
	F:=Filtered(SG,H->Size(H) mod p<>0);;
	SortBy(F,Size);;
	if p=5 then Remove(F);fi; #for p=5 the last semidirect product cannot be handled as a SmallGroup, so one has to do it by hand.
	for H in F do
		G:=SemidirectProduct(H,GF(p)^2);
		if AbelianInvariantsMultiplier(H)<>[] then
			S:=SchurCover(G);
			S:=Image(IsomorphismPermGroup(S));
		else S:=G; fi;
		ct:=CharacterTable(S);
		pb:=PrimeBlocks(ct,p);
		mct:=ct mod p;
		Z:=CenterOfGroupAlgebra(S, p);
		cis:=CentralIdempotentsOfAlgebra(Z);
		for i in [1..Length(pb.defect)] do
			Q:=DecompositionMatrix(mct,i);
			if i>1 then sc:=IdGroup(S);else sc:=IdGroup(G);fi;
			for ci in cis do
				if Dimension(Ideal(Z, [ci]))=Size(Q) then dim:=DimLoewyLayers(Ideal(Z, [ci]));fi;
			od;
			Print(IdGroup(H)," ",IdGroup(G)," ", sc," ",Size(Q)," ",Size(Q[1])," ",dim," ",Length(Positions(ElementaryDivisorsMat(TransposedMat(Q)*Q), 1)),"\n");
		od;
	od;
end;

#Lists all cyclic normal p-blocks in a similar manner as above for D=C_p^n
ListCyclicBlocks:=function(p, n)
	local D, A, I, S, ct, mct, Q, pb, Z, dim;
	D:=CyclicGroup(p^n);
	A:=AutomorphismGroup(D);
	for I in AllSubgroups(A) do
		if Size(I) mod p=0 then continue;fi;
		S:=SemidirectProduct(I, D);
		ct:=CharacterTable(S);
		mct:=ct mod p;
		Q:=DecompositionMatrix(mct);
		pb:=PrimeBlocks(ct, p);
		Z:=CenterOfGroupAlgebra(S, p);
		dim:=DimLoewyLayers(Z);
		Print(IdGroup(I)," ",IdGroup(S)," ",Size(Q)," ",Size(Q[1])," ",pb.defect," ",dim," ",Length(Positions(ElementaryDivisorsMat(TransposedMat(Q)*Q), 1)),"\n");
	od;
end;
