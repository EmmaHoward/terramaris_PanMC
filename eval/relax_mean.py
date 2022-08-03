import iris

paths = ["/gws/nopw/j04/terramaris/emmah/coupled_N1280/kpp_relaxation/2003/",
"/gws/nopw/j04/terramaris/for_tape/kpp_relaxation/2005/",
"/gws/nopw/j04/terramaris/for_tape/kpp_relaxation/2007/",
"/gws/nopw/j04/terramaris/for_tape/kpp_relaxation/2009/",
"/gws/nopw/j04/terramaris/emmah/coupled_N1280/kpp_relaxation/2012/",
"/gws/nopw/j04/terramaris/for_tape/kpp_relaxation/2014/",
"/gws/nopw/j04/terramaris/for_tape/kpp_relaxation/2015/2015",
"/gws/nopw/j04/terramaris/for_tape/kpp_relaxation/2016_6day/2016",
"/gws/nopw/j04/terramaris/for_tape/kpp_relaxation/2017/2017",
"/gws/nopw/j04/terramaris/for_tape/kpp_relaxation/2018/2018"]

cz = iris.Constraint(z=lambda x: -5>=x>=-6)

out = []
for path in paths:
    data=iris.load(path+"*temp*.nc",cz)
    data = iris.cube.CubeList([cube[:-1] for cube in data])
    iris.util.equalise_attributes(data)
    if path[-4:]=='2016':
       data = iris.cube.CubeList([cube for cube in data if cube.shape[0]==6])
    data = data.concatenate_cube()
    out.append(data)


ct = iris.Constraint(t= lambda tt: tt.point.month in [12,1,2])
out=iris.cube.CubeList(out)

