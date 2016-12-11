TARGET = mpi
OBJS += init.o
OBJS += boundary.o
OBJS += uvp.o
OBJS += compute.o
OBJS += main.o

CC = mpicc
#CFLAGS = -Wall -Werror -DDEBUG -g # debug flags
CFLAGS = -Wall -Werror -g -O3 # release flags
CFLAGS += -MMD -MP
LDFLAGS += $(LIBS)

default:	$(TARGET)
all:		$(TARGET)

$(TARGET):	$(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

DEPS = $(OBJS:%.o=%.d)
-include $(DEPS)

clean:
	rm $(TARGET) $(OBJS) $(DEPS) latedays.qsub.* || true
