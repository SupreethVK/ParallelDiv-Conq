all: paradivandconq.out divandconq.out

paradivandconq.out: paralleldivandconq.c
	gcc -o paradivandconq.out paralleldivandconq.c -pthread
	
divandconq.out: divandconq.c
	gcc -o divandconq.out divandconq.c
