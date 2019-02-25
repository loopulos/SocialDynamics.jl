#definition of particles
struct particle
    pos::Array{Float64,1} #postitions
    pgrid::Array{Float64,1} #positions in grid
    vel::Array{Float64,1} #velocities
    mass::Float64 #mass
    indbox::Int64 #the box where the particle is
    tp::Int64
    function particle(dim::Array{Int64,1}, m::Float64, tp::Int64)
        this = new()
        this.pos = rand(2) .* dim
        this.pgrid = this.pos[:]
        this.vel = [rand()*rand([-1,1]),rand()*rand([-1,1])]
        this.mass = m
        this.indbox = ceil(this.pos[1]) + dim[1] * (ceil(this.pos[2])-1)
        this.tp = tp
        return this
    end
    function particle(x::Int64,dim::Array{Int64,1}, m::Float64, tp::Int64)
        this = new()
        this.pos = [rand()+(x-1),rand()* dim[2]]
        this.pgrid = this.pos[:]
        this.vel = [rand()*rand([-1,1]),rand()*rand([-1,1])]
        this.mass = m
        this.indbox = ceil(this.pos[1]) + dim[1] * (ceil(this.pos[2])-1)
        this.tp = tp
        return this
    end
end
#definition of a box
struct box
    vel::Array{Float64,1} # the velocity of the box.
    np::Array{Int64,1}#  not really sure if is neccesary to have the number of particles of each box.
    ind::Int64
    function box(m::Int64, x::Int64)
        this = new()
        this.vel = zeros(2)
        this.np = zeros(Int64,m)
        this.ind = x
        return this
    end
end
