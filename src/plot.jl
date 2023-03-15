function plotDag(g::AbstractGraph,gParams::Dict;verbose=false)
    nodelabels = Vector{String}(undef,0)
    nodestyles = Dict{Int64,String}()
    edgelabels = Dict{Tuple{Int64,Int64},String}()
    edgestyles = Dict{Tuple{Int64,Int64},String}()
    # apply default vertex styles
    for v in 1:nv(g)
        if haskey(props(g,v),:displayName)
            verbose ? println("adding $v") : nothing
            nlab = join([v,props(g,v)[:displayName]],":")
            push!(nodestyles,v=>"fill=white")
        elseif haskey(props(g,v),:intPubTitles)
            verbose ? println("adding $v") : nothing
            nlab = join([v,props(g,v)[:intPubTitles]],":")
            push!(nodestyles,v=>"fill=white")
        elseif haskey(props(g,v),:entId)
            verbose ? println("adding $v") : nothing
            nlab = join([v,props(g,v)[:entId]],":")
            push!(nodestyles,v=>"fill=white")
        else
            verbose ? println("adding $v") : nothing
            nlab = "$v"
            push!(nodestyles,v=>"fill=white")
        end
        if length(nlab) > gParams[:nodelabelmax]
            nlab = string(nlab[1:gParams[:nodelabelmax]],"...")
            push!(nodelabels, nlab)
        else
            push!(nodelabels, nlab)
        end
    end

    if haskey(gParams,:vfilter)
        for s in 1:length(gParams[:vfilter])
            sty = collect(gParams[:vfilter][s])
            verts = filterVertices(g,sty[1][1],sty[1][2][1])
            for v in verts
                nodestyles[v] = sty[1][2][2]
            end
        end
    end

    if haskey(gParams,:highlightnodes)
        for s in gParams[:highlightnodes]
            nodestyles[s] = "fill=green!50"
        end
    end

    # apply default edge styles
    for e in edges(g)
        s = src(e)
        d = dst(e)
        et = (s,d)
        verbose ? println("found edge $s -> $d") : nothing
        push!(edgestyles,et=>"black")
        push!(edgelabels,et=>"")
    end


    # process custom key styles and filter
    if haskey(gParams,:efilter)
        for s in 1:length(gParams[:efilter])
            sty = collect(gParams[:efilter][s])
            edges = filterEdges(g,sty[1][1],sty[1][2][1])
            for e in edges
                edgestyles[(src(e),dst(e))] = sty[1][2][2]
                elab = sty[1][2][3](e,edgelabels)
                if length(elab) > gParams[:edgelabelmax]
                    elab = string(elab[1:gParams[:edgelabelmax]],"...")
                    edgelabels[(src(e),dst(e))] = elab
                else
                    edgelabels[(src(e),dst(e))] = elab
                end
            end
        end
    end

    if haskey(gParams,:highlightedges)
        for e in edges(g)
            s = src(e)
            d = dst(e)
            if any(in.(s,eval(gParams[:highlightedges]))) && any(in.(d,eval(gParams[:highlightedges])))
                et = (s,d)
                verbose ? println("highlighting edge $s -> $d") : nothing
                push!(edgestyles,et=>"red,line width=1em")
                push!(edgelabels,et=>"")
            end
        end
    end

    ## escape "_" which generates a latex compile error because tries to do \text{..., _, ...} instead of \text{..., \_, ...}
    nodelabels = replace.(nodelabels,  r"(\_)" => s"\\\1")
    map!(l->replace(l,  r"(\_)" => s"\\\1"), values(edgelabels))


    verbose ? println("rendering plot") : nothing
    gp = TikzGraphs.plot(DiGraph(g),gParams[:lt],nodelabels,
                    node_style="draw,rounded corners",
                    node_styles=nodestyles,
                    edge_style="black",
                    edge_styles=edgestyles,
                    edge_labels=edgelabels,
                    options=gParams[:opt])
    if split(gParams[:fname],".")[end] == "tex"
        TikzPictures.save(TEX(gParams[:fname]), gp)
    elseif split(gParams[:fname],".")[end] == "svg"
        TikzPictures.save(SVG(gParams[:fname]), gp)
    elseif split(gParams[:fname],".")[end] == "pdf"
        TikzPictures.save(PDF(gParams[:fname]), gp)
    else
        verbose ? println("unknown output format: ",split(gParams[:fname],".")[end]) : nothing
        throw(error())
    end
    gp
end

# create a colorscale for the plot based on max value orthologs across all nodes
function plotDagExp(g::AbstractGraph,gParams::Dict;verbose=false)
    nodequants = Vector{Union{Missing,Float64}}(undef,0)
    nodelabels = Vector{String}(undef,0)
    nodestyles = Dict{Int64,String}()
    edgelabels = Dict{Tuple{Int64,Int64},String}()
    edgestyles = Dict{Tuple{Int64,Int64},String}()

    if haskey(gParams,:colorscheme)
        cscheme = gParams[:colorscheme]
    else
        cscheme = ColorSchemes.plasma
    end


    # get the max value of all expressed orthologs at each node
    for v in 1:nv(g)
        quantmax = missing
        if haskey(props(g,v),:orthoDist)
            if haskey(props(g,v)[:orthoDist],gParams[:orthodist])
                vorth = props(g,v)[:orthoDist][gParams[:orthodist]][:members]
                if length(vorth) > 0
                    verbose ? println("length of vorth is ",length(vorth)) : nothing
                    # for o in 1:length(vorth)
                    #     if haskey(vorth[o],:rnaSeq)
                    #     end
                    # end
                    if length(collect(skipmissing([haskey(vorth[o],:rnaSeq) ? vorth[o][:rnaSeq][:value] : missing for o in 1:length(vorth)]))) > 0
                        quantmax = reduce(max,skipmissing([haskey(vorth[o],:rnaSeq) ? vorth[o][:rnaSeq][:value] : missing for o in 1:length(vorth)]))
                        verbose ? println(quantmax) : nothing
                    else
                        verbose ? println("no expression of ortholog") : nothing
                    end
                else
                    verbose ? println("vorth is empty") : nothing
                end
            else
                verbose ? println("orthoDist=",gParams[:orthodist]," missing") : nothing
            end
        else
            verbose ? println("all orthoDist missing") : nothing
        end
        push!(nodequants,quantmax)
    end
    # derive a log colorscale for the expression values
    if length(collect(skipmissing(nodequants))) == 0
        verbose ? println("no orthologs at dist=",gParams[:orthodist]," detected, exiting") : nothing
        return nothing
    end
    scaledquants = copy(nodequants)
    if gParams[:scale] == :log
        logquants = log.(collect(skipmissing(nodequants)).+1)
        α = 1/maximum(logquants)
        for v in 1:nv(g)
            if !ismissing(scaledquants[v])
                scaledquants[v] = log(scaledquants[v])*α
            end
        end
    elseif gParams[:scale] == :ident
        α = 1/maximum(collect(skipmissing(scaledquants)))
        for v in 1:nv(g)
            if !ismissing(scaledquants[v])
                scaledquants[v] = scaledquants[v]*α
            end
        end
    end


    for v in 1:nv(g)
        # set node styles as expression, when possible
        if !ismissing(scaledquants[v])
            scaled = get(cscheme,scaledquants[v])
            rgbstr = string("fill={rgb:red,",scaled.r,
                                 ";green,",scaled.g,
                                 ";blue,",scaled.b,"}")
            push!(nodestyles,v=>rgbstr)
        else
            push!(nodestyles,v=>"fill=white")
        end

        # set node labels
        if haskey(props(g,v),:displayName)
            verbose ? println("adding $v") : nothing
            nlab = join([v,props(g,v)[:displayName]],":")
        elseif haskey(props(g,v),:intPubTitles)
            verbose ? println("adding $v") : nothing
            nlab = join([v,props(g,v)[:intPubTitles]],":")
        elseif haskey(props(g,v),:entId)
            verbose ? println("adding $v") : nothing
            nlab = join([v,props(g,v)[:entId]],":")
        else
            verbose ? println("adding $v") : nothing
            nlab = "$v"
        end
        if length(nlab) > gParams[:nodelabelmax]
            nlab = string(nlab[1:gParams[:nodelabelmax]],"...")
            push!(nodelabels, nlab)
        else
            push!(nodelabels, nlab)
        end
    end

    if haskey(gParams,:vfilter)
        for s in 1:length(gParams[:vfilter])
            sty = collect(gParams[:vfilter][s])
            verts = filterVertices(g,sty[1][1],sty[1][2][1])
            for v in verts
                nodestyles[v] = sty[1][2][2]
            end
        end
    end

    # apply default edge styles
    for e in edges(g)
        s = src(e)
        d = dst(e)
        et = (s,d)
        verbose ? println("found edge $s -> $d") : nothing
        push!(edgestyles,et=>"black")
        push!(edgelabels,et=>"")
    end

    # process custom key styles and filter
    if haskey(gParams,:efilter)
        for s in 1:length(gParams[:efilter])
            sty = collect(gParams[:efilter][s])
            edges = filterEdges(g,sty[1][1],sty[1][2][1])
            for e in edges
                edgestyles[(src(e),dst(e))] = sty[1][2][2]
                elab = sty[1][2][3](e)
                if length(elab) > gParams[:edgelabelmax]
                    elab = string(elab[1:gParams[:edgelabelmax]],"...")
                    edgelabels[(src(e),dst(e))] = elab
                else
                    edgelabels[(src(e),dst(e))] = elab
                end
            end
        end
    end

    ## escape "_" which generates a latex compile error because tries to do \text{..., _, ...} instead of \text{..., \_, ...}
    nodelabels = replace.(nodelabels,  r"(\_)" => s"\\\1")
    map!(l->replace(l,  r"(\_)" => s"\\\1"), values(edgelabels))

    gp = TikzGraphs.plot(DiGraph(g),gParams[:lt],nodelabels,
                    node_style="draw,rounded corners",
                    node_styles=nodestyles,
                    edge_style="black,line width=0.25em",
                    options=gParams[:opt])

    if split(gParams[:fname],".")[end] == "tex"
        TikzPictures.save(TEX(gParams[:fname]), gp)
    elseif split(gParams[:fname],".")[end] == "svg"
        TikzPictures.save(SVG(gParams[:fname]), gp)
    elseif split(gParams[:fname],".")[end] == "pdf"
        TikzPictures.save(PDF(gParams[:fname]), gp)
    else
        verbose ? println("unknown output format: ",split(gParams[:fname],".")[end]) : nothing
        throw(error())
    end
    gp
end
