import { MainMenu } from './menu.js';
import { OrbitGame } from './orbitGame.js';

export class LevelMenu {
  constructor(canvas, ctx, switchScene) {
    this.canvas = canvas;
    this.ctx = ctx;
    this.switchScene = switchScene;

    const starColors = ['white', '#ffd700', '#add8e6'];
    this.stars = Array.from({ length: 100 }, () => ({
      x: Math.random() * this.canvas.width,
      y: Math.random() * this.canvas.height,
      alpha: Math.random(),
      speed: 0.005 + Math.random() * 0.01,
      drift: (Math.random() - 0.5) * 0.01, 
      color: starColors[Math.floor(Math.random() * starColors.length)],
    }));

    this.shootingStar = null;
    this.shootingTimer = 0;
    this.shootingCooldown = 300 + Math.random() * 300; 
  }

  update() {
    this.stars.forEach(star => {
      star.alpha += star.speed;
      if (star.alpha > 1 || star.alpha < 0) {
        star.speed *= -1;
        star.alpha = Math.max(0, Math.min(1, star.alpha));
      }

      star.x += star.drift;
      if (star.x < 0) star.x = this.canvas.width;
      if (star.x > this.canvas.width) star.x = 0;
    });

    this.shootingTimer++;
    if (!this.shootingStar && this.shootingTimer > this.shootingCooldown) {
      this.spawnShootingStar();
      this.shootingTimer = 0;
      this.shootingCooldown = 300 + Math.random() * 300;
    }

    if (this.shootingStar) {
      this.shootingStar.x += this.shootingStar.dx;
      this.shootingStar.y += this.shootingStar.dy;
      this.shootingStar.alpha -= 0.02;

      if (this.shootingStar.alpha <= 0) {
        this.shootingStar = null;
      }
    }
  }

  spawnShootingStar() {
    const startX = Math.random() * this.canvas.width;
    const startY = Math.random() * this.canvas.height / 2;
    const angle = Math.random() * Math.PI / 4; 
    const speed = 10;

    this.shootingStar = {
      x: startX,
      y: startY,
      dx: speed * Math.cos(angle),
      dy: speed * Math.sin(angle),
      alpha: 1.0,
    };
  }

  render() {
    const ctx = this.ctx;
    ctx.fillStyle = 'black';
    ctx.fillRect(0, 0, this.canvas.width, this.canvas.height);

    this.stars.forEach(star => {
      ctx.beginPath();
      ctx.globalAlpha = star.alpha;
      ctx.fillStyle = star.color;
      ctx.arc(star.x, star.y, 1.5, 0, Math.PI * 2);
      ctx.fill();
    });
    ctx.globalAlpha = 1.0;

    if (this.shootingStar) {
      const s = this.shootingStar;
      ctx.beginPath();
      ctx.strokeStyle = `rgba(255,255,255,${s.alpha})`;
      ctx.lineWidth = 2;
      ctx.moveTo(s.x, s.y);
      ctx.lineTo(s.x - s.dx * 2, s.y - s.dy * 2);
      ctx.stroke();
    }

    ctx.textAlign = 'center';
    ctx.fillStyle = 'white';
    ctx.font = '20px Courier New';
    ctx.fillText('1. Full Mission', this.canvas.width / 2, this.canvas.height / 2 - 20);
    ctx.fillText('2. Final Approach', this.canvas.width / 2, this.canvas.height / 2 + 20);
    ctx.fillText('3. Docking', this.canvas.width / 2, this.canvas.height / 2 + 60);
    ctx.fillStyle = 'yellow';
    ctx.fillText('Press ESC to Quit', this.canvas.width / 2, this.canvas.height / 2 + 100);
  }

  onKeyPress(e) {
    if (['Digit1', 'Numpad1'].includes(e.code)) {
      this.switchScene(new OrbitGame(this.canvas, this.ctx, this.switchScene, 1));
    } else if (['Digit2', 'Numpad2'].includes(e.code)) {
      this.switchScene(new OrbitGame(this.canvas, this.ctx, this.switchScene, 2));
    } else if (['Digit3', 'Numpad3'].includes(e.code)) {
      this.switchScene(new OrbitGame(this.canvas, this.ctx, this.switchScene, 3));
    } else if (e.code === 'Escape') {
      this.switchScene(new MainMenu(this.canvas, this.ctx, this.switchScene));
    }
  }
}
