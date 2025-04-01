import { TutorialView } from './tutorial.js';
import { LevelMenu } from './levelMenu.js';

export class MainMenu {
  constructor(canvas, ctx, switchScene) {
    this.canvas = canvas;
    this.ctx = ctx;
    this.switchScene = switchScene;

    // Spectral palette
    const starColors = ['white', '#ffd700', '#add8e6']; // white, gold, light blue

    this.stars = Array.from({ length: 100 }, () => ({
      x: Math.random() * this.canvas.width,
      y: Math.random() * this.canvas.height,
      alpha: Math.random(),
      speed: 0.005 + Math.random() * 0.01,
      drift: (Math.random() - 0.5) * 0.01, // subtle parallax movement
      color: starColors[Math.floor(Math.random() * starColors.length)],
    }));

    // Shooting star setup
    this.shootingStar = null;
    this.shootingTimer = 0;
    this.shootingCooldown = 300 + Math.random() * 300; // ~5–10 seconds

    this.logo = new Image();
    this.logo.src = 'assets/logo.png';
  }

  update() {
    this.stars.forEach(star => {
      // Twinkle
      star.alpha += star.speed;
      if (star.alpha > 1 || star.alpha < 0) {
        star.speed *= -1;
        star.alpha = Math.max(0, Math.min(1, star.alpha));
      }

      // Drift
      star.x += star.drift;
      if (star.x < 0) star.x = this.canvas.width;
      if (star.x > this.canvas.width) star.x = 0;
    });

    // Shooting star logic
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
    const angle = Math.random() * Math.PI / 4; // shallow angle
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
    const canvas = this.canvas;

    ctx.fillStyle = 'black';
    ctx.fillRect(0, 0, canvas.width, canvas.height);

    // Stars
    this.stars.forEach(star => {
      ctx.beginPath();
      ctx.globalAlpha = star.alpha;
      ctx.fillStyle = star.color;
      ctx.arc(star.x, star.y, 1.5, 0, Math.PI * 2);
      ctx.fill();
    });
    ctx.globalAlpha = 1.0;

    // Shooting star
    if (this.shootingStar) {
      const s = this.shootingStar;
      ctx.beginPath();
      ctx.strokeStyle = `rgba(255,255,255,${s.alpha})`;
      ctx.lineWidth = 2;
      ctx.moveTo(s.x, s.y);
      ctx.lineTo(s.x - s.dx * 2, s.y - s.dy * 2);
      ctx.stroke();
    }

    // Logo
    if (this.logo.complete) {
      const scale = 0.25;
      const w = this.logo.width * scale;
      const h = this.logo.height * scale;
      ctx.drawImage(this.logo, canvas.width / 2 - w / 2, canvas.height / 2 - 250, w, h);
    }

    // Title and UI
    ctx.textAlign = 'center';
    ctx.fillStyle = '#318CE7';
    ctx.font = '40px Courier New';
    ctx.fillText('Lunar Lander Challenge', canvas.width / 2, canvas.height / 2 - 60);

    ctx.font = '20px Courier New';
    ctx.fillStyle = 'white';
    ctx.fillText('Press ENTER to Start', canvas.width / 2, canvas.height / 2);
    ctx.fillStyle = 'yellow';
    ctx.fillText('Press T for Tutorial', canvas.width / 2, canvas.height / 2 + 40);
    ctx.fillStyle = 'white';
    ctx.fillText('Press ESC to Quit', canvas.width / 2, canvas.height / 2 + 80);
  }

  onKeyPress(e) {
    switch (e.code) {
      case 'Enter':
      case 'NumpadEnter':
        this.switchScene(new LevelMenu(this.canvas, this.ctx, this.switchScene));
        break;
      case 'KeyT':
        this.switchScene(new TutorialView(this.canvas, this.ctx, this.switchScene));
        break;
      case 'Escape':
        console.log('Exit (placeholder)');
        break;
    }
  }
}
