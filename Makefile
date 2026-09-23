FC = gfortran
FFLAGS = -O3 
LDFLAGS =

TARGET = test_cg
OBJS = precision_mod.o test_funcs_mod.o cg_mod.o test_cg.o

.PHONY: all
all: $(TARGET)

$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) -o $@ $(LDFLAGS)

precision_mod.o: precision_mod.f90
	$(FC) $(FFLAGS) -c $< -o $@

test_funcs_mod.o: test_funcs_mod.f90 precision_mod.o
	$(FC) $(FFLAGS) -c $< -o $@

cg_mod.o: cg_mod.f90 precision_mod.o
	$(FC) $(FFLAGS) -c $< -o $@

test_cg.o: test_cg.f90 cg_mod.o test_funcs_mod.o precision_mod.o
	$(FC) $(FFLAGS) -c $< -o $@

.PHONY: run
run: $(TARGET)
	./$(TARGET)

.PHONY: clean
clean:
	rm -f *.o *.mod $(TARGET)

.PHONY: cleanall
cleanall: clean
	rm -f result.txt

.PHONY: rebuild
rebuild: clean all
