# 编译输出文件加尾缀 "_VERSION"
VERSION = 0.1
# 源文件 .c 文件的所在的名为 src 的目录
SOURCE = $(wildcard ./src/*.c)
# 增加一个名为 src1 的源文件所在目录
#SOURCE += $(wildcard ./src1/*.c)

# 引用文件 .h 文件的所在的名为 inc 的目录
INCLUEDS = -I./inc
# 增加一个名为 inc1 的引用文件所在目录
#INCLUEDS += -I ./inc1

# 编译输出目标文件的文件名
TARGET    = pd
# 编译输出目标文件所在的文件夹名
OUTPUTDIR = obj
# 编译器全名
CC        = gcc

# 以下都是固定的不用动

CFLAGS  = -Wall -g
OBJECT  = $(patsubst %.c, %.o, $(SOURCE))

$(TARGET): $(OBJECT)
	@mkdir -p $(OUTPUTDIR)/
	$(CC) $^ $(CFLAGES) -o $(OUTPUTDIR)/$(TARGET).exe
	@cp $(OUTPUTDIR)/$(TARGET).exe ./

%.o: %.c
	$(CC) $(INCLUEDS) $(CFLAGES) -c $< -o $@

.PHONY:clean
clean:
	@rm -rf $(OBJECT) $(OUTPUTDIR)/
