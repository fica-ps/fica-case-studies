module fastICA
    export correlation_matrix, covariance_matrix, eigen_decomposition, sing_val_decomp, whiten, fast_ica
    using LinearAlgebra
    using Statistics
    

    correlation_matrix(mat::Array{Float64 ,2})::Array{Float64 ,2} = (mat .- mean(mat,dims=1)) ./ std(mat,dims=1)   
    
    covariance_matrix(corrmat::Array{Float64 ,2})::Array{Float64 ,2} = corrmat |> length |> float |> n -> corrmat' * corrmat / (n  - 1.)

    sing_val_decomp(mat::Array{Float64 ,2})::SVD = mat |> correlation_matrix |> svd
    
    #Returns: a whitened matrix using one of the available whitening techniques(PCA or ZCA)
    #Receives: 
    #mat -> the matrix to whiten 
    #pca -> optional boolean which is used to determine the whitening technique
    function whiten(mat::Array{Float64 ,2},pca::Bool = true)::Array{Float64 ,2}
        U,S,_ = sing_val_decomp(mat)
        values = broadcast(+,  eps(0.3), S)
        pca_w = (Diagonal(values .^ (-1/2)) * U') * mat
        if(pca)
            return pca_w
        end
        #zca whitening
        return U * pca_w
    end 

    #Fast Ica deflation algorithm
    function fast_ica_def(maxiter::Int64, nic::Int64,X::Array{Float64,2},tol::Float64,W::Array{Float64,2},alpha::Float64)::Array{Float64 , 2}     
        #return W
        retW = zeros(size(W))        
        #normalize weight vector to unity i.e. vector sum = 1
        #src : http://www.measurement.sk/2011/Patil.pdf pg 119 
        for i = 1:nic
            wp = W[i,:,]
            wp = normalize!(wp,1)
            W[i,:,] = wp
            #to-do create aux func
            if (i > 1)
                t = zeros(size(wp))
                for u = 1:(i-1)
                   k = sum(wp * retW[u,:,]')
                   t = t + k * retW[u,:,]
                end
                wp = wp - t
            end
            iter = 0
            chg = 0
            converge = false
            while !converge && iter < maxiter
                iter+=1
                wx = wp' * X
                gwx = tanh.(alpha * wx)
                gwx = vcat(gwx,gwx)
                xgwx = X .* gwx
                v1 = vec(mapslices(mean, xgwx, dims = 2))
                gdwx = alpha * (1 .- tanh.(alpha * wx).^2)
                v2 = mean(gdwx) * wp
                w1 = v1 - v2
                #to-do create aux func
                if (i > 1) 
                    t = w1
                    t = zeros(size(w1))
                    for u = 1:(i-1)
                        k = sum(w1 * retW[u,:,]')
                        t = t + k * retW[u,:,]
                    end
                    w1 = w1 - t
                end
                w1 = w1 / sqrt(sum(w1.^2))   
                println("wp = $wp")
                println("w1 = $w1")
                println("W for iter $iter = $W")
                #check for convergence
                #eg: https://github.com/JuliaStats/MultivariateStats.jl/blob/master/src/ica.jl#L97 ln 124 and 125
                chg =  maximum(abs.(abs.(sum(w1 .* wp)) .- 1.0))
                println("Change for iter $iter = $chg")              
                wp = w1
                converge = ( chg < tol )
            end
            retW[i,:,] = wp
        end
        return retW
    end

   #Returns: the resultant ICA model, an instance of the type ICA
    #Receives:
    #maxiter -> number of iterations for the main loop
    #nic     -> number of independent components
    #X       -> the data matrix, must be whitened(use the whiten function)
    #tol     -> tolerable change of the weights at convergence
    #
    #optional :
    #alpha   -> must be between [1,2]
    #whiten  -> whiten the data
    function fast_ica(maxiter::Int64, nic::Int64,X::Array{Float64,2},tol::Float64, alpha::Float64 = 1. , whiten ::Bool=false)::Array{Float64,2}
        #validate arguments
        m,n = size(X)        
        m > 1 || error("There must be more than one samples, n > 1.")
        maxiter > 1 || error("maxiter must be greater than 1.")
        tol > 0 || error("tol must be positive.")
        alpha >= 1 && alpha <= 2|| error("alpha must be in between 1 and 2")
        K = 0
        comp = min(n,m)
        if(nic > comp)
          nic = comp
        end
        X1 = X
        if(whiten)
            X = (X .- mean(X,dims=1))'
            V = X * X'/n
            s = svd(V)
            D = Diagonal(1 ./ sqrt.(s.S))
            K = D * s.U'
            X1 = K * X
        end
        # initialize weights of size n with random values
        #src : https://en.wikipedia.org/wiki/FastICA 
        W = rand(comp,comp)
        a = fast_ica_def(maxiter,nic,X1,tol,W,alpha)
        println("K=$K")
        w = a * K
        return w * X
    end
    

#=
        a = fast_ica_def(maxiter,nic,X1,tol)
        w = a * K
        return w * X
    ### PARALLEL FAST_ICA IMPLEMENTATION USING LOGCOSH
    
    function fast_ica(maxiter::Int64, nic::Int64,X::Array{Float64,2},tol::Float64, alpha::Int64 = 1)::Array{Float64 , 2}
        #validate arguments
        p,n = size(X)
        n > 1 || error("There must be more than one samples, n > 1.")
        maxiter > 1 || error("maxiter must be greater than 1.")
        tol > 0 || error("tol must be positive.")
        alpha >= 1 && alpha <= 2|| error("alpha must be in between 1 and 2")
        
        if(p > n)
            X = X'
        end
        if(nic > min(n,p))
            nic = min(n,p)
        end
        
        # initialize weights of size n with random values
        #src : https://en.wikipedia.org/wiki/FastICA 
        W = rand(nic,nic)
        #old W
        oldW = zeros(nic,nic)        
        #normalize weight vector to unity i.e. vector sum = 1
        #src : http://www.measurement.sk/2011/Patil.pdf pg 119 
        for i = 1:nic
            wp = view(W, i,:,)
            wp = normalize!(wp,1)
        end
        
        sW = svd(W)
        W = sW.U * Diagonal(1 ./ sW.S) * sW.U' * W
        #main loop
        # src: 
        #http://www.measurement.sk/2011/Patil.pdf pg 119 Fixed Point algorithm for ICA
        #https://www.cs.helsinki.fi/u/ahyvarin/papers/NN00new.pdf pg 15
        t = 0
        chg = 0
        converge = false
        while !converge && t < maxiter
             t+=1
            #store previous iteration W
            copyto!(oldW,W)
            
            #based on CRAN implementation ( parallel using logcosh aprox. to neg-entropy function )
            #src : https://cran.r-project.org/web/packages/fastICA/fastICA.pdf
            wx = W * X
            gwx = tanh.(alpha * wx)
            v1 = gwx * X'/p
            gdwx = alpha * (1 .- gwx.^2)
            v2 = Diagonal(vec(mapslices(mean, gdwx, dims = 2))) * W
            W = v1 - v2
            sW1 = svd(W)
            W = sW1.U * Diagonal( 1 ./ sW.S ) * sW1.U' * W    
            println("W for iter $t = $W")
            #check for convergence
            #eg: https://github.com/JuliaStats/MultivariateStats.jl/blob/master/src/ica.jl#L97 ln 124 and 125
            chg =  maximum(abs.(abs.(diag(W*oldW')) .- 1.0))
            println("Change for iter $t = $chg")
            converge = ( chg < tol )
        end
        converge || throw("Convergence was not possible with tol = $tol, maxiter= $maxiter, last change = $chg")
        return W
    end
    
    ######
    function fast_ica(maxiter::Int64, nic::Int64,X::Array{Float64,2},tol::Float64)::Array{Float64 , 2}
        #validate arguments
        n,m = size(X)
        n > 1 || error("There must be more than one samples, n > 1.")
        maxiter > 1 || error("maxiter must be greater than 1.")
        tol > 0 || error("tol must be positive.")
        chg = 0
        # initialize weights of size n with random values
        #src : https://en.wikipedia.org/wiki/FastICA 
        W = rand(n,nic)
        #old W
        oldW = Array{Float64}(undef,n,nic)
        for p = 1:nic
            #normalize weight vector to unity i.e. vector sum = 1
            #src : http://www.measurement.sk/2011/Patil.pdf pg 119
            wp = view(W, :, p)
            wp = normalize!(wp,1)
            iter = 0
            converge = false
            while  iter < maxiter
                iter+=1
                #store previous iteration W
                copyto!(oldW,W)
                t = X'*wp;
                g = t.^3; 
                dg = 3*t.^2; 
                wp = X*g/m-mean(dg)*wp;
                
                #symmetric decorrelation
                #src : https://www.cs.helsinki.fi/u/ahyvarin/papers/NN00new.pdf pg 15
                #eg : https://github.com/JuliaStats/MultivariateStats.jl/blob/master/src/ica.jl ln 120
                copyto!(W, W * (W'W).^(-1/2))     
                #check for convergence
                #eg: https://github.com/JuliaStats/MultivariateStats.jl/blob/master/src/ica.jl#L97 ln 124 and 125
                chg = 1 - maximum(abs.(abs.(dot.(W,oldW)) .- 1.0))
                
                if( chg <= tol )
                    return W
                end
            end
        end
        throw("Convergence was not possible with tol = $tol, maxiter= $maxiter \n last iteration change = $chg")
    end

 =#
end