classdef Motivation < HybridSystem
    
    properties(SetAccess = immutable)
        theta
        s
        gammac
        gammad
    end

    methods 
        function this = Motivation(parameters)
            state_dim = 11; % (z,thetahat,tau_s,tau,q)
            this = this@HybridSystem(state_dim);
            
            this.theta = parameters.theta;
            this.s = parameters.s;
            this.gammac = parameters.gammac;
            this.gammad = parameters.gammad;
        end

        function xdot = flowMap(this, x, t, j)
            psi = x(1:2);
            thetahat = x(3:4);
            tau_s = x(5);
            tau = x(6);
            q = x(7);
            thetahatc = x(8:9);
            thetahatd = x(10:11);
            
            psidot = [cos(t); 0];
            thetahatdot = 0*thetahat;
            tau_sdot = 1;
            taudot = 1;
            qdot = 0;
            thetahatcdot = 0*thetahatc;
            thetahatddot = 0*thetahatd;

            xdot = [psidot; thetahatdot; tau_sdot; taudot; qdot; thetahatcdot; thetahatddot];
        end

        function xplus = jumpMap(this, x, t, j)
            psi = x(1:2);
            thetahat = x(3:4);
            tau_s = x(5);
            tau = x(6);
            q = x(7);
            thetahatc = x(8:9);
            thetahatd = x(10:11);
            
            inD1 = tau_s >= this.s;     % take a sample during flows
            inD2 = tau >= 2*pi;         % jump due to jump set
            inD3 = q == 1;              % jump due to jump set
            
            if inD1        % take a sample during flows
                psiplus = psi;
                tau_splus = 0;
                tauplus = tau;
                qplus = q;
                
                y = psi'*this.theta;
                thetahatplus = thetahat + this.s*this.gammac*psi*(y - psi'*thetahat);
                
                thetahatcplus = thetahatc + this.s*this.gammac*psi*(y - psi'*thetahatc);
                thetahatdplus = thetahatd;
                
            elseif inD2 || inD3    % jump due to the jump set
                if inD2
                    psiplus = [0.5; 1];
                    qplus = 1;
                    tauplus = 0;
                else
                    psiplus = [0; 0];
                    qplus = 0;
                    tauplus = tau;
                end
                tau_splus = tau_s;
                
                yplus = psiplus'*this.theta;
                thetahatplus = thetahat + this.gammad*psiplus*(yplus - psiplus'*thetahat)/(1 + this.gammad*norm(psiplus)^2);
                
                thetahatcplus = thetahatc;
                thetahatdplus = thetahatd + this.gammad*psiplus*(yplus - psiplus'*thetahatd)/(1 + this.gammad*norm(psiplus)^2);
                
            else
                zplus = z;
                thetahatplus = thetahat;
                tauplus = tau;
                thetahatplusc = thetahatc;
                thetahatplusd = thetahatd;
            end
            
            xplus = [psiplus; thetahatplus; tau_splus; tauplus; qplus; thetahatcplus; thetahatdplus];
        end
        
        function inC = flowSetIndicator(this, x)
            inC = 1;
        end

        function inD = jumpSetIndicator(this, x)
            psi = x(1:2);
            thetahat = x(3:4);
            tau_s = x(5);
            tau = x(6);
            q = x(7);
            thetahatc = x(8:9);
            thetahatd = x(10:11);
            
            inD1 = tau_s >= this.s;     % take a sample during flows
            inD2 = tau >= 2*pi;         % jump due to jump set
            inD3 = q == 1;              % jump due to jump set
            inD = inD1 || inD2 || inD3;
        end
    end
end
