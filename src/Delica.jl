module Delica

using Unitful, Tessen, Statistics, LinearAlgebra

export polycontour, sectorcontour, box

"""
```julia
polycontour(vertices,[fillet])
```
Create a `Tessen.Contour` representing a polygon with `vertices` in the xy plane.
Optional fillet parameter can be used to round the corners. All arguments should have
units of `Unitful.Length`
"""
function polycontour(vertices::Vector{<:Vector{<:Unitful.Length}},fillet::Unitful.Length)
    @assert all(length.(vertices) .== 2) "All vertices should consist of two coordinates"
    #create a vector of edge endpoints
    vertcopy = copy(vertices)
    push!(vertcopy,vertices[1])
    edges = map(1:length(vertices)) do i
        [vertcopy[i],vertcopy[i+1]]
    end
    #if iszero(fillet) then we're done
    if iszero(fillet)
        return Contour([LineEdge(points...) for points in edges])
    end
    #do some bookkeeping
    edgecopy = copy(edges)
    push!(edgecopy,edges[1])
    vertedges = map(1:length(edges)) do i
        (leftedge=edgecopy[i],rightedge=edgecopy[i+1])
    end
    #now for each entry in vertedges we need to calculate:
    #leftpoint, the start of the fillet,
    #fillet, the ArcEdge representing the fillet
    #rightpoint, the end of the fillet
    #the points can then be connected with lineedges to make the contour
    fillets = map(vertedges) do ve
        @assert ve.leftedge[2] == ve.rightedge[1]
        #get the address of 'our' vertex
        vcoords = ve.leftedge[2]
        #get vectors parallel to leftedge and rightedge
        vectors = [ve.leftedge[1] - vcoords, ve.rightedge[2] - vcoords]
        #get the angle between each of these vectors and the x axis
        angles = map(vectors) do v
            atan(reverse(v)...)
        end
        #the center of the fillet lies on a line which bisects the angle formed
        #by our edges
        cangle = mean(angles)
        #get half of the angle of this vertex
        vangle = abs(cangle-angles[1])
        #get the length of the line segment connecting 'our' vertex to the center
        #of the fillet
        cdist = fillet / sin(vangle)
        #this is all the information we need to get the coordinates of the fillet center
        cdisp = cdist*[cos(cangle),sin(cangle)]
        #we always want to fillet concave angles, so the dot product of cdisp with either of
        #our `vectors` should be positive
        dot(cdisp,vectors[1])
        if dot(cdisp,vectors[1]) < (zero(cdisp[1])^2)
            cdisp *= -1
            vangle = pi - vangle
        end
        c = vcoords + cdisp
        #we can also get the distance between the vertex and the tangency point with the fillet
        tdist = cdist * cos(vangle)
        #this is enough to calculate the locations of leftpoint and rightpoint
        points = map(angles) do a
            vcoords + tdist*[cos(a),sin(a)]
        end
        #we can now get the angle between the x axis and a vector connecting c to each of
        #`points`
        pangles = map(points) do p
            thisangle = atan(reverse(p - c)...)
            #make sure thisangle is in the range [-2pi,2pi]
            if abs(thisangle) > 2pi
                thisangle %= 2pi
            end
            #now make it between [0 and 2pi]
            (thisangle >= 0) ? thisangle : thisangle + 2pi
        end
        #we want to draw the shorter of the two possible arcs connecting `points`
        anglediff = pangles[2] - pangles[1]
        #if pangles is in the correct order, anglediff will be less than pi and positive
        #or greater than pi and negative
        correctorder = if anglediff > 0
            abs(anglediff) <= pi
        else
            abs(anglediff) >= pi
        end
        if !correctorder
            reverse!(pangles)
        end
        fillet
        (leftpoint = points[1], fillet = ArcEdge(c,fillet,pangles...), rightpoint = points[2])
    end
    filletcopy = copy(fillets)
    push!(filletcopy,fillets[1])
    alledges = vcat(map(1:length(fillets)) do i
                        fillet
                        [filletcopy[i].fillet,
                         LineEdge(filletcopy[i].rightpoint,filletcopy[i+1].leftpoint)]
                    end...)
    Contour(alledges)
end

polycontour(vertices) = polycontour(vertices, 0u"μm")

"""
```julia
annularcontour(c,r1,r2 [,startangle,stopangle; rounded=false])
```
Create a closed contour representing a sector of an annulus centered on point `c`, with inner
and outer radii `r1` and `r2`. If `startangle` and `stopangle` are not provided, the contour
will be a complete annulus. If the optional keyword argument `rounded` is set to `true`, the
endcaps of the contour will be semicircles rather than straight lines.
"""
function sectorcontour end

#if there's no startangle and stopangle this is dead simple
function sectorcontour(c::Vector{<:Unitful.Length},r1::Unitful.Length,r2::Unitful.Length)
    map([r1,r2]) do r
        ArcEdge(c,r,0,2pi)
    end |> Contour
end

#otherwise...
function sectorcontour(c::Vector{<:Unitful.Length},r1::Unitful.Length,r2::Unitful.Length,
                       startangle::Number, stopangle::Number;
                       rounded=false)
    #making our primary arcs are easy
    mainarcs = map([r1,r2]) do r
        ArcEdge(c,r,startangle,stopangle)
    end
    #caps are lines or arcs depending on the value of rounded
    caps = if rounded
        #arcs
        #get the radius of the midline
        rmid = mean([r1,r2])
        #thickness of the annulus
        tann = r2 - r1
        #the caps need to be swept in opposite directions, we will use `forward` to keep track
        map(zip([startangle,stopangle],[false,true])) do (a,forward)
            #get the center of the cap
            ccap = c + rmid*[cos(a), sin(a)]
            angles = [a,a+pi]
            if !forward
                reverse!(angles)
            end
            ArcEdge(ccap,tann/2,angles...)
        end
    else
        #lines
        map([startangle,stopangle]) do a
            cappoints = map([r1,r2]) do r
                #get an endpoint of the endcap line
                c + r*[cos(a), sin(a)]
            end
            LineEdge(cappoints...)
        end
    end
    #make our closed contour
    Contour(vcat(mainarcs,caps))
end

"""
```julia
box(length, width, height, dslice;
    [chamfer,fillet])
```
Build a 3D box with the provided dimensions. `dslice` is the maximum
allowable slicing distance (the number of required layers and lines
will be rounded up). The dimension corresponding to `length` will be
along the x axis. `chamfer` should be a matrix with size (2,2)
specifying the angle with which the edges of the box should be tapered.
chamfer[1,:] specifies tapering normal to the 'length' dimension, chamfer[2,:]
specifies tapering normal to `width`. If `chamfer` is provided the box has a
crossection of ``length × width`` at the z coordinate corresponding to
the center of the object. If `fillet` is provided, each crossection is filleted.
This function will build a box centered on `[0,0,0]`. `Tessen.rotate` and
`Tessen.translate` can be used to make new objects with the same dimensions in
different locations.
"""
function box(length::Unitful.Length, width::Unitful.Length,
             height::Unitful.Length,dslice::Unitful.Length;
             chamfer=zeros(Float64,2,2),fillet = 0u"µm")
    @assert size(chamfer) == (2,2)
    #get the z position of all of our layers
    numz = ceil(Int,height/dslice) #minimum value for dslice
    zpos=range(-height/2,
               height/2, length = numz) |> collect
    #build all of the slices
    slices = map(1:Base.length(zpos)) do i
        #apply the chamfer
        zoffset=zpos[i]
        #amount of material added to each edge
        #negative sign so positive chamfers correspond to shapes
        #without overhang
        added = -1*(zoffset * tan.(chamfer))
        #the change in the center of the crosssection
        center = (added[:,2] - added[:,1]) / 2
        (lprime,wprime) = [length,width] + sum(added,dims=2)
        #calculate our corner positions
        rawcorners = [[-lprime/2, -wprime/2],
                      [lprime/2,  -wprime/2],
                      [lprime/2,  wprime/2],
                      [-lprime/2, wprime/2]]
        corners = [rc + center for rc in rawcorners]
        thisslice = Slice([polycontour(corners,fillet)])
        (zpos[i] => thisslice)
    end
    Block(slices...)
end

end # module Delica
