using ThreadingUtilities

struct RGFWrap{F, Args, T}
    f::F
end
@inline args_type(f::RGFWrap{F,Args}) where {F, Args} = Args

function (f::RGFWrap{F, Args, T})(p::Ptr{UInt}) where {Args, F, T}
    _, (argsptr, resptr) = ThreadingUtilities.load(p, Tuple{Ptr{Args}, Ptr{T}}, 2*sizeof(UInt))
    res = f.f(unsafe_load(argsptr)...)
    unsafe_store!(resptr, res)
    nothing
end

@inline function fptr(f::RGFWrap)
    @cfunction($f, Cvoid, (Ptr{UInt},))
end

struct Funcall{F, Args, outT}
    f::RGFWrap{F, Args, outT}
    args::Args
end

(f::Funcall)() = f.f.f(f.args...)

output_type(f::Funcall{F, Args, outT}) where {F, Args, outT} = outT

@inline function Funcall(f::F, args::Args, outT=Base.return_types(f, Args)[1]) where {F,Args}
    Funcall{F, Args, outT}(RGFWrap{F, Args, outT}(f), args)
end

function setup_call!(p, f, args, resref)
    fp = Base.unsafe_convert(Ptr{Cvoid}, fptr(f))
    resptr = Base.unsafe_convert(Ptr{eltype(resref)}, resref)
    argptr = Base.unsafe_convert(Ptr{args_type(f)}, args)
    offset = ThreadingUtilities.store!(p, fp, sizeof(UInt))
    ThreadingUtilities.store!(p, (argptr, resptr), offset)
end

@inline function launch_call!(tid, f::RGFWrap, args, resref)
    p = ThreadingUtilities.taskpointer(tid)
    while true
        if ThreadingUtilities._atomic_cas_cmp!(p, ThreadingUtilities.SPIN, ThreadingUtilities.STUP)
            setup_call!(p, f, args, resref)
            @assert ThreadingUtilities._atomic_cas_cmp!(p, ThreadingUtilities.STUP, ThreadingUtilities.TASK)
            return
        elseif ThreadingUtilities._atomic_cas_cmp!(p, ThreadingUtilities.WAIT, ThreadingUtilities.STUP)
            setup_call!(p, f, args, resref)
            @assert ThreadingUtilities._atomic_cas_cmp!(p, ThreadingUtilities.STUP, ThreadingUtilities.LOCK)
            ThreadingUtilities.wake_thread!(tid % UInt)
            return
        end
        ThreadingUtilities.pause()
    end
end
@generated function spawn_fetch(fs::NTuple{N,Any}, ::Val{nt}, g=tuple) where {nt, N}

    ts = 1:nt
    batches = map(Iterators.partition(1:N, nt)) do batch
         # one batch of spawns
         launches = map(ts) do t
             if t > length(batch)
                 return
             end
             i = batch[t]
             if t == nt
                 :(rs[$i][] = fs[$i]())
             else
                 :(launch_call!($(ts[t]),
                                fs[$i].f,
                                argrefs[$i],
                                rs[$i]))
             end
         end

         waits = map(ts) do t
             if t > length(batch)
                 return
             end
             i = batch[t]
             if t == nt
                 :(rs[$i][] = fs[$i]())
             else
                 :(ThreadingUtilities.__wait($(ts[t])))
             end
         end

         quote
             $(launches...)
             $(waits...)
         end
     end

    quote
        rs = Base.@ntuple $N i->Ref{output_type(fs[i])}()
        argrefs = Base.@ntuple $N i->Ref(fs[i].args)
        GC.@preserve rs argrefs begin
            $(batches...)
        end
        return Base.@ntuple $N i->rs[i][]
    end
end

function spawn_fetch_nonstatic(fs, g=tuple)
    rs = map(f->Ref{output_type(f)}(), fs)
    argrefs = map(f->Ref(f.args), fs)

    nt = Threads.nthreads()
    ts = 1:nt-1
    GC.@preserve rs argrefs begin
        for batch in Iterators.partition(1:length(fs),
                                         length(ts))
            for (t, i) in enumerate(batch)
                launch_call!(ts[t], fs[i].f, argrefs[i], rs[i])
            end
            for (t, i) in enumerate(batch)
                ThreadingUtilities.__wait(t)
            end
        end
    end
    return g(map(r->r[], rs)...)
end


const DEBUGBUF = zeros(UInt, 64)
@generated function spawn_fetch_debug(fs::NTuple{N,Any}, ::Val{nt}, g=tuple) where {nt, N}

    ts = 1:nt
    batches = map(Iterators.partition(1:N, nt)) do batch
         # one batch of spawns
         launches = map(ts) do t
             if t > length(batch)
                 return
             end
             i = batch[t]
             if t == nt
                 :(rs[$i][] = fs[$i]())
             else
                 quote
                    setup_call!(p,
                                fs[$i].f,
                                argrefs[$i],
                                rs[$i])
                    fs[$i].f(p)
                end
             end
         end

         quote
             $(launches...)
         end
     end

    quote
        rs = Base.@ntuple $N i->Ref{output_type(fs[i])}()
        argrefs = Base.@ntuple $N i->Ref(fs[i].args)
        p = Base.unsafe_convert(Ptr{UInt}, DEBUGBUF)
        GC.@preserve rs argrefs begin
            $(batches...)
        end
        return Base.@ntuple $N i->rs[i][]
    end
end

#=
function spawn_fetch_debug(fs, g=tuple)
    rs = map(f->Ref{output_type(f)}(), fs)
    argrefs = map(f->Ref(f.args), fs)

    nt = Threads.nthreads()
    ts = 1:nt-1
    arr = zeros(UInt, 64)
    GC.@preserve arr fs rs argrefs begin
        for batch in Iterators.partition(1:length(fs),
                                         length(ts))
            for (t, i) in enumerate(batch)
                fill!(arr, 0)
                p = Base.unsafe_convert(Ptr{UInt}, arr)
                setup_call!(p, fs[i].f, argrefs[i], rs[i])
                fs[i].f(p)
            end
        end
    end
    return g(map(r->r[], rs)...)
end
=#
