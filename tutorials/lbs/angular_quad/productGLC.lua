--[[ @doc
# Product Quadrature of Gauss-Legendre-Chebyshev type.
___
## Create the Angular Quadrature
  + ```nazimu``` is the number of subdivisions in **one*** quadrant of the the equatorial plane,
  + ```npolar``` is the number of **positive** polar cosines.

Hence, there will be $nazimu \times npolar$ directions per octant.

In 2D XY geometry, we only use 4 octants instead of 8 octants.
-- @end ]]
--################################################
nazimu = 2
npolar = 4
pquad = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,nazimu,npolar)
chiOptimizeAngularQuadratureForPolarSymmetry(pquad, 4.0*math.pi)

--[[ @doc
### Here, we retrieve the quadrature data in lua to print it
-- @end]]
qdata = chiGetProductQuadrature(pquad)

print("")
print(string.format("+---------+-------------+------------+------------+"))
print(string.format("|direction|    weight   | cos(theta) | phi (deg.) |"))
print(string.format("+---------+-------------+------------+------------+"))
for d,v in pairs(qdata) do
  s = string.format("| %5d   |    %6.4f   |   %7.4f  |   %7.3f  |",d, v.weight, math.cos(v.polar), v.azimuthal*180.0/math.pi)
  print(s)
end
print(string.format("+---------+-------------+------------+------------+"))

--[[ @doc
### Printing quadrature data to a file for subsequent plotting.

You can plot the directions of the quadrature using the [Python script plot_ang_quad.py](./plot_ang_quad.py) found in this same folder.

A sample plot is shown below:
![Directions](images/ang_quad_plot.png)

-- @end]]
-- Opens a file in write mode
filename = "qdata_na" .. tostring(nazimu) .. "_np" .. tostring(npolar) .. ".csv"
fhandle = io.open (filename , "w")

-- sets the default output file
io.output(fhandle)

-- write to the file
for d,v in pairs(qdata) do
  io.write(v.weight, ", ", math.cos(v.polar), ", ", v.azimuthal, "\n")
end

-- closes the open file
io.close(fhandle)

--[[ @doc
### Another manner to print out data
The lua code below is current unused, see commented line.
-- @end ]]
--################################################
-- Print contents of `tbl`, with indentation.
-- `indent` sets the initial level of indentation.
-- Usage: tprint(qdata)
function tprint (tbl, indent)
  if not indent then indent = 0 end
  for k, v in pairs(tbl) do
    formatting = string.rep("  ", indent) .. k .. ": "
    if type(v) == "table" then
      print(formatting)
      tprint(v, indent+1)
    elseif type(v) == 'boolean' then
      print(formatting .. tostring(v))
    else
      print(formatting .. v)
    end
  end
end
