SWI= -c -O
airfoil: airfoil.o airfoil_flow.o
	f77 -o airfoil airfoil.o airfoil_flow.o
airfoil.o: airfoil.f airfoil_com.f
	f77 $(SWI) airfoil.f
airfoil_flow.o: airfoil_flow.f airfoil_com.f
	f77 $(SWI) airfoil_flow.f


