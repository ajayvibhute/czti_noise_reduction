#this make file is for AstroSat CZTI Pipeline
#!/bin/csh
SRC_PATH=src

all:
	make -C $(SRC_PATH)/cztscience2event
	make -C $(SRC_PATH)/cztpha2energy
	make -C $(SRC_PATH)/cztbunchclean
	make -C $(SRC_PATH)/cztgaas
	make -C $(SRC_PATH)/cztgtigen
	make -C $(SRC_PATH)/cztdatasel
	make -C $(SRC_PATH)/cztpixclean	
	make -C $(SRC_PATH)/cztevtclean
	make -C $(SRC_PATH)/cztflagbadpix
	make -C $(SRC_PATH)/cztdpigen
	make -C $(SRC_PATH)/cztimage
	make -C $(SRC_PATH)/cztbindata
	make -C $(SRC_PATH)/cztrspgen

	@echo "----------------------------------------"
	@echo "Compilation of All modules of CZTI pipeline completed successfully"
	@echo "----------------------------------------"




clean:
	
	make -C $(SRC_PATH)/cztscience2event clean
	make -C $(SRC_PATH)/cztpha2energy clean
	make -C $(SRC_PATH)/cztbunchclean clean
	make -C $(SRC_PATH)/cztgaas clean
	make -C $(SRC_PATH)/cztgtigen clean
	make -C $(SRC_PATH)/cztdatasel clean
	make -C $(SRC_PATH)/cztpixclean clean
	make -C $(SRC_PATH)/cztevtclean clean
	make -C $(SRC_PATH)/cztflagbadpix clean
	make -C $(SRC_PATH)/cztdpigen clean
	make -C $(SRC_PATH)/cztimage clean
	make -C $(SRC_PATH)/cztbindata clean
	make -C $(SRC_PATH)/cztrspgen clean 


distclean:
	rm  bin/*
	
		


