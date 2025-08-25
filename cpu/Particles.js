class Particles {

    constructor() {
        this._particles = [];
    }

    Add(p) {
        this._particles.push(p);
    }

    AddRange(particleArray) {
        if (!Array.isArray(particleArray)) {
            throw new Error("AddRange expects an array of particles.");
        }

        for (let p of particleArray) {
            this.Add(p);
        }
    }

    Delete(p) {
        this._particles = this._particles.filter(item => item.id !== p.id);
    }

    DeleteRange(particleArray) {
        if (!Array.isArray(particleArray)) {
            throw new Error("DeleteRange expects an array of particles.");
        }

        for (let p of particleArray) {
            this.Delete(p);
        }
    }

    GetAll() {
        return this._particles;
    }

    Clear() {
        this._particles.length = 0;
    }

    Count() {
        return this._particles.length;
    }

    FindById(id) {
        return this._particles.find(p => p.id === id);
    }

    Draw(ctx, transformer){
        this._particles.forEach(p => p.Draw(ctx, transformer));
    }
}
