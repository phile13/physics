class Particle {
    static ID = 0;
    static TWO_PI = Math.PI * 2;
    

    constructor(type, x, y, vx, vy){
        this.id = Particle.ID++;
        this.type = this.SetType(type);
        this.radius = 1e-16;
        if(this.type === "Photon"){
            this.radius = .5e-17;
        }

        this.current = {x : x, y : y, vx : vx || 0, vy : vy || 0};
        this.next = {x : x, y : y, vx : vx || 0, vy : vy || 0};
        this.fx = 0;
        this.fy = 0;
        this.frame_data = {x : this.next.x , y : this.next.y , q : this.charge_unitless, color : this.fill_color, type : this.type, vx : this.current.vx , vy  : this.current.vy };

        this.age = this.life_span * Math.random();
        this.radiation_buildup = 0;

        const v2 = vx * vx + vy * vy;
        if(v2 >= PhysicsEngine.C_SQR || this.type === "Photon"){
            this.px = 0;
            this.py = 0;
        }
        else{
            const gamma_times_mass = (1 / Math.sqrt(1 - Math.min(v2 / (PhysicsEngine.C_SQR),1.00000))) * this.mass;
            // Relativistic momentum
            this.px = gamma_times_mass * vx;
            this.py = gamma_times_mass * vy;
        }
    }

    Draw(ctx, transformer){
        this.frame_data = {x : this.next.x , y : this.next.y , q : this.charge_unitless , color : this.fill_color, type : this.type, vx : this.next.vx || 0 , vy : this.next.vy || 0};
        let [pX,pY,pR] = transformer(this.next.x, this.next.y, this.radius);
        ctx.beginPath();
        ctx.arc(pX,pY,pR, 0, Particle.TWO_PI);
        ctx.fillStyle = this.fill_color;
        ctx.fill();
    }

    Clear(){
        [this.next,this.current] = [this.current, this.next];
        this.fx = this.fy = 0;
    }

    SetType(type){
        [   
            this.is_force,
            this.mass,
            this.charge,
            this.life_span,
            this.fill_color, 
            this.anti_particle
        ] 
        = Particle.PARTICLE_TYPES[type];
        
        this.charge_unitless = this.charge / 1.602176634e-19;
        this.mass_sqr = this.mass * this.mass;
        this.mass_sqr_times_c_sqr = this.mass_sqr * PhysicsEngine.C_SQR;
        return type;
    }

    

    static PARTICLE_TYPES = {
        Electron :[false,9.109383e-31,-1.602176634e-19,6.6e28,"blue","Positron"],
        Positron :[false,9.109383e-31,1.602176634e-19,6.6e28,"orange","Electron"],
        Muon :[false,1.883531e-28,-1.602176634e-19,2.196e-22,"pink","AntiMuon"],
        AntiMuon :[false,1.883531e-28,1.602176634e-19,2.196e-6,"palegreen","Muon"],
        Tau :[false,3.16754e-27,-1.602176634e-19,2.903e-13,"red","AntiTau"],
        AntiTau :[false,3.16754e-27,1.602176634e-19,2.903e-13,"green","Tau"],
        Photon : [true,0,0,6.6e28,"black",""],
        VirtualPhoton : [true,0,0,8.5e-24,"rgba(0, 68, 255, 0.4)"

        ]
    };

    static CheckForOverlap(p1, p2) {
        // Extract initial and final positions
        const x1_0 = p1.current.x, y1_0 = p1.current.y;
        const x1_1 = p1.next.x,    y1_1 = p1.next.y;
        const x2_0 = p2.current.x, y2_0 = p2.current.y;
        const x2_1 = p2.next.x,    y2_1 = p2.next.y;
    
        const r1 = p1.radius || 0.1, r2 = p2.radius || 0.1;
        const R = r1 + r2;
    
        // Relative motion: treat p2 as stationary and p1 moving with relative velocity
        const dx0 = x1_0 - x2_0;
        const dy0 = y1_0 - y2_0;
        const dx1 = x1_1 - x2_1;
        const dy1 = y1_1 - y2_1;
    
        const dvx = dx1 - dx0;
        const dvy = dy1 - dy0;
    
        // solve for |(dx0,dy0) + t*(dvx,dvy)|^2 = R^2
        const a = dvx*dvx + dvy*dvy;
        const b = 2 * (dx0*dvx + dy0*dvy);
        const c = dx0*dx0 + dy0*dy0 - R*R;
    
        let overlap = false;
        let t = null;
    
        if (a === 0) {
            // No relative movement, overlap only if already intersecting
            overlap = c <= 0;
            t = overlap ? 0 : null;
        } else {
            const disc = b*b - 4*a*c;
            if (disc >= 0) {
                const sqrt_disc = Math.sqrt(disc);
                const t1 = (-b - sqrt_disc) / (2*a);
                const t2 = (-b + sqrt_disc) / (2*a);
                // Select the earliest valid collision time within [0, 1]
                let candidates = [t1, t2].filter(time => time >= 0 && time <= 1);
                if (candidates.length > 0) {
                    overlap = true;
                    t = Math.min(...candidates);
                }
            }
        }
    
        if (overlap) {
            const x1 = p1.current.x + (p1.next.x - p1.current.x) * t;
            const y1 = p1.current.y + (p1.next.y - p1.current.y) * t;
            const x2 = p2.current.x + (p2.next.x - p2.current.x) * t;
            const y2 = p2.current.y + (p2.next.y - p2.current.y) * t;
    
            const nx = x1 - x2, ny = y1 - y2;
            const nrm = Math.hypot(nx, ny) || 1;
            const collisionNormal = {x: nx / nrm, y: ny / nrm};
    
            return {
                timeOfImpact: t,
                contactPoint1: {x: x1, y: y1},
                contactPoint2: {x: x2, y: y2},
                collisionNormal: collisionNormal,
                penetrationDepth: Math.max(0, (p1.radius + p2.radius) - nrm)
            };
        }
    
        return null;
    }
    
}
