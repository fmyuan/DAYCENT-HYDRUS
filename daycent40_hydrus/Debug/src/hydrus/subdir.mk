################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
F_SRCS += \
../src/hydrus/HYDRUS_DayCent.f \
../src/hydrus/INPUT_DayCent.f \
../src/hydrus/MATERIAL.f \
../src/hydrus/OUTPUT_DAYCENT.f \
../src/hydrus/SETBC_DAYCENT.f \
../src/hydrus/SINK.f \
../src/hydrus/TIME.f \
../src/hydrus/WATFLOW.f 

OBJS += \
./src/hydrus/HYDRUS_DayCent.o \
./src/hydrus/INPUT_DayCent.o \
./src/hydrus/MATERIAL.o \
./src/hydrus/OUTPUT_DAYCENT.o \
./src/hydrus/SETBC_DAYCENT.o \
./src/hydrus/SINK.o \
./src/hydrus/TIME.o \
./src/hydrus/WATFLOW.o 


# Each subdirectory must supply rules for building sources it contributes
src/hydrus/%.o: ../src/hydrus/%.f
	@echo 'Building file: $<'
	@echo 'Invoking: GNU Fortran Compiler'
	$(FC) -fno-underscoring -O0 -g -Wall -c -fmessage-length=0 -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

src/hydrus/HYDRUS_DayCent.o: ../src/hydrus/HYDRUS_DayCent.f

src/hydrus/INPUT_DayCent.o: ../src/hydrus/INPUT_DayCent.f

src/hydrus/MATERIAL.o: ../src/hydrus/MATERIAL.f

src/hydrus/OUTPUT_DAYCENT.o: ../src/hydrus/OUTPUT_DAYCENT.f

src/hydrus/SETBC_DAYCENT.o: ../src/hydrus/SETBC_DAYCENT.f

src/hydrus/SINK.o: ../src/hydrus/SINK.f

src/hydrus/TIME.o: ../src/hydrus/TIME.f

src/hydrus/WATFLOW.o: ../src/hydrus/WATFLOW.f


