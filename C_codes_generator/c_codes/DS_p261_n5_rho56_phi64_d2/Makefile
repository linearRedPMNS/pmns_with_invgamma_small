OPTIONS= -O3 -Wall -lgmp
#~ OPTIONS= -O3 -Wall -g 

CC= gcc

OPTIONS+=$(MACROS)


all : main

main : main.o pmns_useful_functs.o pmns_arith_ops.o
	$(CC) -o main main.o pmns_useful_functs.o pmns_arith_ops.o $(OPTIONS)


main.o : main.c
	$(CC) -c main.c $(OPTIONS)


pmns_useful_functs.o : pmns_useful_functs.c 
	$(CC) -c pmns_useful_functs.c $(OPTIONS)


pmns_arith_ops.o : pmns_arith_ops.c 
	$(CC) -c pmns_arith_ops.c $(OPTIONS)


clean :
	rm -f *.o *~
	rm -f main

