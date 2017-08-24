k:=FiniteField(9);
G:=SmallGroup(72,40);
B:=BasicAlgebraOfGroupAlgebra(G,k);
quiver, rels, crels:=QuiverAndRelations(B);
print "Quiver of the principal block:";
print quiver;
print "Quiver relations:";
print crels;

print "\nRadical layers of projective indecomposable modules";
for i:=1 to 5 do
	i;
	IsomorphismTypesOfRadicalLayers(ProjectiveModule(B,i));
end for;
