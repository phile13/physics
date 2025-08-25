class Simulation {
    static MinSize = 1e-17;
    static NumPixelsLightCoversInAFrame = 1;
    static DistanceLightCoversInAFrame = Simulation.NumPixelsLightCoversInAFrame * Simulation.MinSize;
    static ParticleForceRingSpacing = 5;
    static ParticleForcesHistorySize = 10000;

    constructor(canvas) {
        this.canvas = canvas;
        this.ctx = canvas.getContext("2d");
        this.particles = new Particles();
        this.physics = new PhysicsEngine((d,c)=>{this.OnParticleEvent(d,c);});
        this.dt = Simulation.DistanceLightCoversInAFrame / 299792458;//light moves 10 pixels per frame, or 1e-17 meters
        this.ulx = 0;
        this.uly = 0;
        this.width = 1000 * Simulation.MinSize;
        this.height = this.width;
        this.tx = this.canvas.width / this.width;
        this.ty = this.canvas.height / this.height;
    }

    Setup(user_defined_particles = []){
        this.particles.AddRange(user_defined_particles);

        let step = 100;
        for(let i = 100; i < 1000; i+=step){
            for(let j = 100; j < 1000; j+=step){
                this.particles.Add(new Particle((Math.random() < .5) ? "Positron" : "Electron", i * Simulation.MinSize, j * Simulation.MinSize, (Math.random()*2-1) * PhysicsEngine.C, (Math.random()*2-1) * PhysicsEngine.C));
            }
        }
        /*//this.particles.Add(new Particle("Electron", 300 * Simulation.MinSize, 100 * Simulation.MinSize, 0, -2.9e8));
        this.particles.Add(new Particle("Electron", 300 * Simulation.MinSize, 200 * Simulation.MinSize, 0, -2.9e8));
        this.particles.Add(new Particle("Electron", 300 * Simulation.MinSize, 300 * Simulation.MinSize, 0, -2.9e8));
        this.particles.Add(new Particle("Electron", 300 * Simulation.MinSize, 400 * Simulation.MinSize, 0, -2.9e8));
        this.particles.Add(new Particle("Electron", 300 * Simulation.MinSize, 500 * Simulation.MinSize, 0, -2.9e8));
        this.particles.Add(new Particle("Electron", 300 * Simulation.MinSize, 600 * Simulation.MinSize, 0, -2.9e8));
        this.particles.Add(new Particle("Electron", 300 * Simulation.MinSize, 700 * Simulation.MinSize, 0, -2.9e8));
        this.particles.Add(new Particle("Electron", 300 * Simulation.MinSize, 800 * Simulation.MinSize, 0, -2.9e8));
        this.particles.Add(new Particle("Electron", 300 * Simulation.MinSize, 900 * Simulation.MinSize, 0, -2.9e8));

        //this.particles.Add(new Particle("Positron", 350 * Simulation.MinSize, 100 * Simulation.MinSize, 0, 0));
        this.particles.Add(new Particle("Positron", 350 * Simulation.MinSize, 200 * Simulation.MinSize, -1e8, 0));
       this.particles.Add(new Particle("Positron",  350 * Simulation.MinSize, 300 * Simulation.MinSize, 1e8, 0));
        this.particles.Add(new Particle("Positron", 350 * Simulation.MinSize, 400 * Simulation.MinSize, -1e8, 0));
        this.particles.Add(new Particle("Positron", 350 * Simulation.MinSize, 500 * Simulation.MinSize, 1e8, 0));
        this.particles.Add(new Particle("Positron", 350 * Simulation.MinSize, 600 * Simulation.MinSize, 1e8, 0));
        this.particles.Add(new Particle("Positron", 350 * Simulation.MinSize, 700 * Simulation.MinSize, -1e8, 0));
        this.particles.Add(new Particle("Positron", 350 * Simulation.MinSize, 800 * Simulation.MinSize, 1e8, 0));
        this.particles.Add(new Particle("Positron", 350 * Simulation.MinSize, 900 * Simulation.MinSize, -1e8, 0));

        this.particles.Add(new Particle("Electron", 400 * Simulation.MinSize, 200 * Simulation.MinSize, 0, 2.9e8));
        this.particles.Add(new Particle("Electron", 400 * Simulation.MinSize, 300 * Simulation.MinSize, 0, 2.9e8));
        this.particles.Add(new Particle("Electron", 400 * Simulation.MinSize, 400 * Simulation.MinSize, 0, 2.9e8));
        this.particles.Add(new Particle("Electron", 400 * Simulation.MinSize, 500 * Simulation.MinSize, 0, 2.9e8));
        this.particles.Add(new Particle("Electron", 400 * Simulation.MinSize, 600 * Simulation.MinSize, 0, 2.9e8));
        this.particles.Add(new Particle("Electron", 400 * Simulation.MinSize, 700 * Simulation.MinSize, 0, 2.9e8));
        this.particles.Add(new Particle("Electron", 400 * Simulation.MinSize, 800 * Simulation.MinSize, 0, 2.9e8));
        this.particles.Add(new Particle("Electron", 400 * Simulation.MinSize, 900 * Simulation.MinSize, 0, 2.9e8));


        //this.particles.Add(new Particle("Positron", 500 * Simulation.MinSize, 100 * Simulation.MinSize, 0, 0));
        //this.particles.Add(new Particle("Electron", 500 * Simulation.MinSize, 500 * Simulation.MinSize, 0, 0));*/
        
        let Ps = this.particles.GetAll();
        this.physics.Init(Ps);
    }

    Run() {
        let Ps = this.particles.GetAll();
        this.physics.Update(Ps, this.dt);
        
        this.ctx.clearRect(0, 0, this.canvas.width, this.canvas.height);

        /*this.physics.DrawForcesHistory(this.ctx, (x, y, radius)=>{
            return [(x-this.ulx)*this.tx, (y-this.uly)*this.ty, radius * this.tx];
        });

        this.physics.DrawTrajectories(this.ctx, (x, y, radius)=>{
            return [(x-this.ulx)*this.tx, (y-this.uly)*this.ty, radius * this.tx];
        });*/

        this.particles.Draw(this.ctx, (x, y, radius)=>{
            return [(x-this.ulx)*this.tx, (y-this.uly)*this.ty, radius * this.tx];
        });
        
        requestAnimationFrame(() => this.Run());
    }

    OnParticleEvent(destroying, creating){
        this.particles.DeleteRange(destroying);
        this.particles.AddRange(creating);
        if(destroying && destroying.length > 0){
            this.physics.CleanupAfterDeletion(destroying);
        }
    }
    
}
