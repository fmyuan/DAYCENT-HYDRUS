################################################################################
# Automatically-generated file. Do not edit!
################################################################################

-include ../makefile.init

RM := rm -rf
FC := /usr/local/gcc-4.8.1/bin/gfortran
CC := /usr/local/gcc-4.8.1/bin/gcc

# All of the sources participating in the build are defined here
-include sources.mk
-include src/hydrus/subdir.mk
-include src/daycent/subdir.mk
-include subdir.mk
-include objects.mk

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
endif

-include ../makefile.defs

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: daycent40_hydrus

# Tool invocations
daycent40_hydrus: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: MacOS X Fortran Linker'
	$(FC)  -o "daycent40_hydrus" $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(OBJS)$(C_DEPS)$(EXECUTABLES) daycent40_hydrus
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:

-include ../makefile.targets
