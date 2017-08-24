LoewyLayers:=function(M, irr)
	LL:=[];
	while Dimension(M) ne 0 do
		layer:=[];
		rad:=JacobsonRadical(M);
		cfs:=CompositionFactors(quo<M | rad>);
		for cf in cfs do
			for i:=1 to #irr do
				if Dimension(AHom(cf, irr[i])) ne 0 then 
					layer:=Append(layer, i-1);
				end if;			
			end for;
		end for;
		Sort(~layer);
		LL:=Append(LL, layer);
		M:=rad;
	end while;
	return LL;
end function;

SocleLayers:=function(M, irr)
	SL:=[];
	SS:=SocleSeries(M);
	for j:=#SS to 1 by -1 do
		layer:=[];
		if j ne 1 then 
			cfs:=CompositionFactors(quo<SS[j] | SS[j-1]>);
		else
			cfs:=CompositionFactors(SS[j]);
		end if;
		for cf in cfs do
			for i:=1 to #irr do
				if Dimension(AHom(cf, irr[i])) ne 0 then
					layer:=Append(layer, i-1);
				end if;
			end for;
		end for;
		Sort(~layer);
		SL:=Append(SL, layer);
	end for;
	return SL;
end function;

G:=Group("A6");
k:=FiniteField(9);
D:=SylowSubgroup(G, 3);
N:=Normalizer(G, D);
IRR:=IrreducibleModules(G, k);
irr:=IrreducibleModules(N, k);
proj:=ProjectiveIndecomposables(N, k);
print "Simple modules of the group";
IRR;
S:=[1,2,3,4];
print "We know the dimensions of the simple modules in the principal block, so we can pick them out as ",S;

print "Loewy structure of the projective indecomposables of the normalizer";
for i:=1 to #irr do
	i-1;
	LoewyLayers(proj[i], irr);
end for;

print "Restricting the simples and calculating their radical and socle layers";
for s in S do
	M:=Restriction(IRR[s], N);
	dec:=Decomposition(M);
	printf "Dim %o Summands %o \n", Dimension(M), #dec;
	for summand in Decomposition(M) do
		Dimension(summand);
		LoewyLayers(summand, irr);
		SocleLayers(summand, irr);
	end for;
end for;
