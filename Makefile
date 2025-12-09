MODULE_TOPDIR = ../..

PGM = i.hyper.rgb

include $(MODULE_TOPDIR)/include/Make/Script.make
include $(MODULE_TOPDIR)/include/Make/Html.make

default: script html
