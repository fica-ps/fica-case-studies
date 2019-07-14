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
        sigma = cov(mat)        
        vals,vecs = eigen(sigma)
        #Sort eigvals and eigvecs by DESCENDING order
        vecs = reverse(vecs, dims=2)
        vals = reverse(vals, dims=1)
        pca_w =  Diagonal(1 ./ sqrt.(vals)) * vecs'
        #U,S,_ = svd(sigma)
        #pca_w =  Diagonal(1 ./ sqrt.(S)) * U'
        if(pca)
            return pca_w
        end
        #zca whitening
        return vecs * pca_w
    end 

    #auxiliary function for logcosh contrast function
    function contrast_func(alpha::Float64,wx)::Tuple{Array{Float64},Array{Float64}}
       gd= tanh.(alpha * wx)
       gdd= alpha * (1 .- tanh.(alpha * wx).^2)
       return (gd,gdd)
    end

    #Fast Ica deflation algorithm
    function fast_ica_def(maxiter::Int64, nic::Int64,X::Array{Float64,2},tol::Float64,W::Array{Float64,2},alpha::Float64)::Array{Float64 , 2}     
        #return W
        retW = zeros(size(W))        
        #src : http://www.measurement.sk/2011/Patil.pdf pg 119 
        for i = 1:nic
            println("Component n $i")
            wp = W[i,:,]
            #to-do create aux func
            if (i > 1)
                t = zeros(size(wp))
                for u = 1:(i-1)
                   k = sum(wp .* retW[u,:,])
                   t = t + k * retW[u,:,]
                end
                wp = wp - t
            end
            normalize!(wp)
            iter = 0
            chg = 0
            converge = false
            while !converge && iter < maxiter
                iter+=1
                wx = wp' * X
                gdx,gddx = contrast_func(alpha,wx)            
                xgdx = X .* gdx
                v1 = vec(mapslices(mean, xgdx, dims = 2))
                v2 = mean(gddx) * wp
                w1 = v1 - v2
                #to-do create aux func
                if (i > 1) 
                    t = zeros(size(w1))
                    for u = 1:(i-1)
                        k = sum(w1 .* retW[u,:,])
                        t = t + k * retW[u,:,]
                    end
                    w1 = w1 - t
                end
                normalize!(w1)
                #check for convergence
                chg =  abs.(abs.(sum(w1 .* wp)) .- 1.0)
                println("Tolerance change for iter $iter = $chg")              
                wp = w1
                converge = ( chg < tol )
            end
            retW[i,:,] = wp
        end
        println("retW = $retW")
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
    function fast_ica(maxiter::Int64, nic::Int64,X::Array{Float64,2},tol::Float64, alpha::Float64 = 1. )::Array{Float64,2}
        #validate arguments
        m,n = size(X)        
        m > 1 || error("There must be more than one samples, n > 1.")
        maxiter > 1 || error("maxiter must be greater than 1.")
        tol > 0 || error("tol must be positive.")
        alpha >= 1 && alpha <= 2|| error("alpha must be in between 1 and 2")
        comp = min(n,m)
        if(nic > comp)
          nic = comp
        end
        # initialize weights of size n with random values
        #src : https://en.wikipedia.org/wiki/FastICA 
        W = randn(comp,comp)
        
        println("Random matrix = $W")
        return fast_ica_def(maxiter,nic,X,tol,W,alpha)
        
    end
    
end