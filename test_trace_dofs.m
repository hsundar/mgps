nelems = [16, 16, 16];
order = 5;

num_elems  = prod(mesh.nelems);
nnf = mesh.refel.nnf;
      ndof = nnf * ( 3*num_elems + mesh.nelems(1)*mesh.nelems(2) + mesh.nelems(2)*mesh.nelems(3) + mesh.nelems(1)*mesh.nelems(3));

      num_bdy = prod(num_elems.*order + 1) - prod(num_elems)*(order-1)^3; 