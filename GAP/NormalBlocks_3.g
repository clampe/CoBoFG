Read("ListNormalBlocks.g");

#Print the blocks we obtained with their inertial quotients, group, schur covers, k(B) and l(B)
p:=3;
Print("List for p=",p,": \n");
Print("Blocks with defect group C_p \n");
ListCyclicBlocks(p, 1);
Print("\nBlocks with defect group C_p^2 \n");
ListCyclicBlocks(p, 2);
Print("\nBlocks with defect group C_p x C_p\n");
ListNormalBlocks(p);
Print("\n");
QUIT_GAP();
