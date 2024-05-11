//! The simplest possible example that does something.
#![allow(clippy::unnecessary_wraps)]

use ggez::{Context, event, GameResult, glam::*, graphics::{self}, mint};
use ggez::conf::{WindowMode, WindowSetup};
use ggez::graphics::{Canvas, DrawMode, FillOptions};
use rand::Rng;

struct MainState {
    grid: Grid,
    ray_origin: Vec2,
    mouse_pos: Vec2,
    ray_dist: f32,
}

impl MainState {
    fn mouse_motion_event(&mut self, _ctx: &mut Context, x: f32, y: f32, _dx: f32, _dy: f32) {
        println!("OKKK");
        self.mouse_pos = Vec2::new(x, y);
    }
    fn new(dim: usize) -> GameResult<MainState> {
        let mut grid: Grid = Grid::new(dim);

        Ok(
            Self {
                grid,
                ray_origin: Vec2::new(0., 0.),
                mouse_pos: Vec2::new(0., 0.),
                ray_dist: 0.,
            }
        )
    }
}

pub const WIN_H: f32 = 1080.;
pub const WIN_W: f32 = 1080.;

pub const cell_w: f32 = WIN_W / (16 as f32);
pub const cell_h: f32 = WIN_H / (16 as f32);

impl event::EventHandler<ggez::GameError> for MainState {
    fn update(&mut self, ctx: &mut Context) -> GameResult {
        let position = ctx.mouse.position();
        self.mouse_pos = Vec2::new(position.x, position.y);
        self.ray_dist = (Vec2::new(0.,0.) - self.mouse_pos).length();
        Ok(())
    }

    fn draw(&mut self, ctx: &mut Context) -> GameResult {
        let mut canvas =
            graphics::Canvas::from_frame(ctx, graphics::Color::from([0.1, 0.2, 0.3, 1.0]));


        for x in 0..self.grid.max_bound.x as usize {
            let line = graphics::Mesh::new_line(
                ctx,
                &[
                    Vec2::new(x as f32 * cell_w, 0.),
                    Vec2::new(x as f32 * cell_w, WIN_H),
                ],
                1.0,
                graphics::Color::from_rgb(255, 255, 0),
            )?;
            canvas.draw(&line, graphics::DrawParam::new());
        }

        for y in 0..self.grid.max_bound.y as usize {
            let line = graphics::Mesh::new_line(
                ctx,
                &[Vec2::new(0., y as f32 * cell_h),
                    Vec2::new(WIN_H, y as f32 * cell_w)
                ],
                1.0,
                graphics::Color::from_rgb(255, 255, 0),
            )?;
            canvas.draw(&line, graphics::DrawParam::new())
        }
        let mouse_line = graphics::Mesh::new_line(
            ctx,
            &[Vec2::new(0.0, 0.0), self.mouse_pos],
            2.0,
            graphics::Color::from_rgb(0, 255, 255),
        )?;

        let dir = (self.mouse_pos - Vec2::new(0., 0.)).normalize();

        let ray = Ray {
            origin: Vec2::new(0., 0.),
            dir,
        };

        let (t0, t1) = (
            0.,
            (self.grid.max_bound - self.grid.min_bound).length()
        );

        traverse(ctx, self, &mut canvas, &ray, t0, f32::clamp(t1, 0., 4.));

        canvas.draw(&mouse_line, graphics::DrawParam::new());
        canvas.finish(ctx)?;

        Ok(())
    }
}

fn traverse(
    ctx: &mut Context,
    state: &mut MainState,
    canvas: &mut Canvas,
    ray: &Ray,
    t0: f32,
    t1: f32,
) {
    let (mut tmin, mut tmax) = (0., 0.);
    let ray_intersects_grid = ray_box_intersect(ray, &state.grid, &mut tmin, &mut tmax, t0, t1);
    if (!ray_intersects_grid) { return; }

    println!("{}", t1);
    println!("{}", tmax);
    (tmin, tmax) = (f32::max(tmin, t0), f32::max(tmax, t1));

    let ray_start: Vec2 = ray.origin + ray.dir * tmin;
    let ray_end: Vec2 = ray.origin + ray.dir * tmax;

    let mut x_idx: i32 = f32::max(1., f32::ceil(ray_start.x - state.grid.min_bound.x)) as i32;
    let x_bound: i32 = f32::max(1., f32::ceil(ray_end.x - state.grid.min_bound.x)) as i32;

    let (step_x, tdx, mut tmax_x) = if ray.dir.x > 0.0 {
        (
            1,
            1. / ray.dir.x,
            tmin + (state.grid.min_bound.x + (x_idx as f32) * 1. - ray_start.x) / ray.dir.x
        )
    } else if ray.dir.x < 0.0 {
        let prev_x_idx = x_idx - 1;
        (
            -1,
            1. / -ray.dir.x,
            tmin + (state.grid.min_bound.x + (prev_x_idx as f32) * 1. - ray_start.x) / ray.dir.x
        )
    } else {
        (
            0,
            tmax,
            tmax
        )
    };

    let mut y_idx: i32 = f32::max(1., f32::ceil(ray_start.y - state.grid.min_bound.y / 1.)) as i32;
    let y_bound: i32 = f32::max(1., f32::ceil(ray_end.y - state.grid.min_bound.y / 1.)) as i32;

    let (step_y, tdy, mut tmax_y): (i32, f32, f32) = if ray.dir.y > 0. {
        (
            1,
            1. / ray.dir.y,
            tmin + (state.grid.min_bound.y + (y_idx as f32) * 1. - ray_start.y) / ray.dir.y
        )
    } else if ray.dir.y < 0.0 {
        let prev_y_idx = y_idx - 1;
        (
            -1,
            1. / -ray.dir.y,
            tmin + (state.grid.min_bound.y + (prev_y_idx as f32) * 1. - ray_start.y) / ray.dir.y
        )
    } else {
        (
            0,
            tmax,
            tmin
        )
    };

    while (x_idx <= 16 && y_idx <= 16) {
        let point = graphics::Mesh::new_circle(
            ctx,
            DrawMode::Fill(FillOptions::DEFAULT),
            mint::Point2 { x: (x_idx as f32) * cell_h, y: (y_idx as f32) * cell_w },
            5.0,
            0.1,
            graphics::Color::from_rgb(255, 0, 0)
        ).unwrap();

        canvas.draw(&point, graphics::DrawParam::new());

        if (tmax_x < tmax_y) {
            x_idx += step_x;
            tmax_x += tdx;
        } else {
            y_idx += step_y;
            tmax_y += tdy;
        }
    }
}


fn ray_box_intersect(
    ray: &Ray,
    grid: &Grid,
    t_min: &mut f32,
    t_max: &mut f32,
    t0: f32,
    t1: f32,
) -> bool {
    let x_inv_dir = 1. / ray.dir.x;
    (*t_min, *t_max) = if (x_inv_dir >= 0.) {
        (
            (grid.min_bound.x - ray.origin.x) * x_inv_dir,
            (grid.max_bound.x - ray.origin.x) * x_inv_dir
        )
    } else {
        (
            (grid.max_bound.x - ray.origin.x * x_inv_dir),
            (grid.min_bound.x - ray.origin.x * x_inv_dir)
        )
    };

    let y_inv_dir = 1. / ray.dir.y;
    let (mut ty_min, mut ty_max) = if y_inv_dir >= 0.0 {
        (
            (grid.min_bound.y - ray.origin.y) * y_inv_dir,
            (grid.max_bound.y - ray.origin.y) * y_inv_dir
        )
    } else {
        (
            (grid.max_bound.y - ray.origin.y) * y_inv_dir,
            (grid.min_bound.y - ray.origin.y) * y_inv_dir
        )
    };

    if *t_min > ty_max || ty_min > ty_max { return false; }
    if ty_min > *t_min { *t_min = ty_min; }
    if ty_max < *t_max { *t_max = ty_max; }

    *t_min < t1 && *t_max > t0
}

struct Grid {
    data: Vec<Vec<u32>>,
    min_bound: Vec2,
    max_bound: Vec2,
}

impl Grid {
    fn new(dim: usize) -> Self {
        let mut data = vec![vec![0; dim]; dim];
        let mut rng = rand::thread_rng();
        for _ in 0..dim {
            let x: usize = rng.gen_range(0..dim) as usize;
            let y: usize = rng.gen_range(0..dim) as usize;
            data[x][y] = 1;
        }
        Self {
            data,
            min_bound: Vec2::new(0., 0.),
            max_bound: Vec2::new(dim as f32, dim as f32),
        }
    }
}

struct Ray {
    origin: Vec2,
    dir: Vec2,
}

pub fn main() -> GameResult {
    let cb = ggez::ContextBuilder::new("super_simple", "ggez")
        .window_setup(WindowSetup::default().title("FUCK"))
        .window_mode(WindowMode::default().dimensions(WIN_W, WIN_H));
    let (mut ctx, event_loop) = cb.build()?;
    let state = MainState::new(16)?;

    event::run(ctx, event_loop, state)
}
