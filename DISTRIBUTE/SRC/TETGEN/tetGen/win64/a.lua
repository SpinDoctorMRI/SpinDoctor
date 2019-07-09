-- Lua script.
p=tetview:new()
p:load_plc("C:/Users/Jing Rebecca LI/Dropbox/SpinDoctor/CODE_JRL/GUI+SRC/voxel.poly")
rnd=glvCreate(0, 0, 500, 500, "TetView")
p:plot(rnd)
glvWait()
