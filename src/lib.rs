#[derive(Debug, Clone, PartialEq, Default)]
pub struct Circle {
    pub x: f32,
    pub y: f32,
    pub radius: f32,
}
/*
math.acos(
    (distance_nodes(cm,cn)**2+(cm.radius+radius)**2-(cn.radius+radius)**2) /
    (2*(cm.radius+radius)*distance_nodes(cm,cn)))
*/
impl Circle {
    pub fn distance_to_origin(&self) -> f32 {
        (self.x * self.x + self.y * self.y).sqrt()
    }
    fn distance(&self, other: &Circle) -> f32 {
        let x = self.x - other.x;
        let y = self.y - other.y;
        (x * x + y * y).sqrt()
    }
    fn intersect(&self, other: &Circle) -> bool {
        // Epsilon Error for floating point errors
        self.distance(other) + 1e-4 < self.radius + other.radius
    }
    fn new_node(radius: f32, c_m: &Circle, c_n: &Circle) -> (Circle, Circle) {
        let dist = c_n.distance(c_m);
        // Non-right angle trig to figure out the angles
        let phi1 = (c_n.y - c_m.y).atan2(c_n.x - c_m.x);
        let phi2 = ((dist.powi(2) + (c_m.radius + radius).powi(2) - (c_n.radius + radius).powi(2))
            / (2.0 * (c_m.radius + radius) * dist))
            .acos();
        // The sin & cos of the angles give the x & y coords on unit circle. Expand and translate.
        let solution_1 = Circle {
            x: c_m.x + (c_m.radius + radius) * (phi1 + phi2).cos(),
            y: c_m.y + (c_m.radius + radius) * (phi1 + phi2).sin(),
            radius,
        };
        // Other solution
        let solution_2 = Circle {
            x: c_m.x + (c_m.radius + radius) * (phi1 - phi2).cos(),
            y: c_m.y + (c_m.radius + radius) * (phi1 - phi2).sin(),
            radius,
        };
        (solution_1, solution_2)
    }

    fn initialise(radii: &[f32]) -> Vec<Circle> {
        match radii.len() {
            0 => Vec::new(),
            1 => vec![Circle {
                x: 0.0,
                y: 0.0,
                radius: radii[0],
            }],
            // Place the two right next to each other
            2 => {
                let r0 = radii[0];
                let r1 = radii[1];
                let ray = r0 + r1;
                vec![
                    Circle {
                        x: ray / 2.0,
                        y: 0.0,
                        radius: r0,
                    },
                    Circle {
                        x: -ray / 2.0,
                        y: 0.0,
                        radius: r1,
                    },
                ]
            }
            // There are 3 or more nodes. We ignore the rest and place the first 3 so they touch.
            _ => {
                let r0 = radii[0];
                let r1 = radii[1];
                let r2 = radii[2];
                let mut front = vec![
                    Circle {
                        x: 0.0,
                        y: 0.0,
                        radius: r0,
                    },
                    Circle {
                        x: r0 + r1,
                        y: 0.0,
                        radius: r1,
                    },
                ];
                {
                    let (sol_1, sol_2) = Self::new_node(r2, &front[0], &front[1]);
                    if sol_1.y < 0.0 {
                        front.push(sol_1);
                    } else {
                        front.push(sol_2)
                    }
                }
                // Center of mass calculation, we put this at the origin
                let (a, b, c) = (r0 + r1, r0 + r2, r1 + r2);
                let s = (a + b + c) / 2.0;
                let delta = (s * (s - a) * (s - b) * (s - c)).sqrt();
                let (l1, l2, l3) = (
                    a + delta / (s - a),
                    b + delta / (s - b),
                    c + delta / (s - c),
                );
                let l_sum = l1 + l2 + l3;
                let (l1, l2, l3) = (l1 / l_sum, l2 / l_sum, l3 / l_sum);
                let x_correction = l1 * front[2].x + l2 * front[1].x + l3 * front[0].x;
                let y_correction = l1 * front[2].y + l2 * front[1].y + l3 * front[0].y;
                front
                    .drain(..)
                    .map(|mut c| {
                        c.x = c.x - x_correction;
                        c.y = c.y - y_correction;
                        c
                    })
                    .collect()
            }
        }
    }
}


#[derive(Debug, Clone, Default)]
pub struct CirclePacker {
    radii: Vec<f32>,
    circles: Vec<Circle>,
    front: Vec<usize>,
}

impl CirclePacker {
    pub fn push(&mut self, radius: f32) {
        self.radii.push(radius);

        if self.radii.len() <= 3 {
            self.circles = Circle::initialise(&self.radii);
            self.front = (0..self.circles.len()).collect();
        } else {
            let mut index = self.front_closest();
            let mut next_index = self.front_get_next_index(index);
            let mut candidate = self.candidate(radius, index, next_index);
            while let Some(i) = self.front_get_intersection(index, &candidate) {
                let (front_dist, back_dist) = self.front_distance_to_node(index, i);

                if front_dist < back_dist {
                    let (break_start, break_end) = self.front_remove(index, i);
                    candidate = self.candidate(
                        radius,
                        break_start,
                        break_end,
                    );
                    index = break_start;
                    next_index = break_end;
                } else {
                    let (break_start, break_end) = self.front_remove(i, next_index);
                    candidate = self.candidate(
                        radius,
                        break_start,
                        break_end,
                    );
                    index = break_start;
                    next_index = break_end;
                }
            }
            self.insert(candidate, next_index);
        }
    }

    fn candidate(&self, radius: f32, index: usize, next_index: usize) -> Circle {
        let (solution_1, solution_2) = Circle::new_node(radius, &self.circles[self.front[index]], &self.circles[self.front[next_index]]);
        let c1 = &self.circles[self.front[index]]; 
        let c2 = &self.circles[self.front[next_index]];
        // The determiniant tells us if a change of basis will flip the orientation or keep it steady. 
        // We want the chosen solution to be on the right side of the vector from c1 to c2's centers. 
        // This means we want the determinant to be negative. 
        let determinate = (solution_1.x - c1.x)*(c2.y - c1.y) - (solution_1.y - c1.y)*(c2.x - c1.x);
        if determinate > 0.0 {
            solution_2
        } else {
            solution_1
        }
    }

    pub fn circles(&self, embedding_circle: Option<&Circle>) -> Vec<Circle> {
        if let (Some(embedding_circle), Some(radius)) = (embedding_circle,self.circles.iter().map(|c| c.distance_to_origin() + c.radius).max_by(|a,b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))) {
            let radius_ratio = embedding_circle.radius / radius;
            self.circles.iter().map(|c| {
                Circle {
                    x: c.x*radius_ratio + embedding_circle.x,
                    y: c.y*radius_ratio + embedding_circle.y,
                    radius: c.radius * radius_ratio,
                }
            }).collect()
        } else {
            self.circles.clone()
        }
    }

    /// Inserts a node
    fn insert(&mut self, node: Circle, front_break_end: usize) {
        let front_index = self.circles.len();
        self.circles.push(node);
        if 0 != front_break_end {
            self.front.insert(front_break_end, front_index);
        } else {
            self.front.push(front_index);
        }
    }

    /// The closest node in the front to the origin
    fn front_closest(&self) -> usize {
        self.front.iter()
                .enumerate()
                .min_by(|(_, m), (_, n)| {
                    self.circles[**m].distance_to_origin()
                        .partial_cmp(&self.circles[**n].distance_to_origin())
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
                .unwrap()
                .0
    }

    /// Wrapping 
    fn front_get_next_index(&self, start_index: usize) -> usize {
        (start_index + 1) % self.front.len()
    }
    /// Find the closest intersecting node in the front to the start index.
    fn front_get_intersection(&self, start_index: usize, candidate: &Circle) -> Option<usize> {
        let mut intersecting_nodes = self.front
            .iter()
            .enumerate()
            .filter_map(|(i, n)| {
                if self.circles[*n].intersect(candidate) {
                    Some(i)
                } else {
                    None
                }
            })
            .collect::<Vec<usize>>();
        intersecting_nodes.sort_by(|a, b| {
            let (a_b, a_f) = self.front_distance_to_node(start_index, *a);
            let a = a_b.min(a_f);
            let (b_b, b_f) = self.front_distance_to_node(start_index, *b);
            let b = b_b.min(b_f);
            a.partial_cmp(&b).unwrap_or(std::cmp::Ordering::Equal)
        });
        intersecting_nodes.first().map(|f| *f)
    }
    /// Distance along the chain from the start to the end.
    fn front_distance_to_node(&self, start_index: usize, end_index: usize) -> (f32, f32) {
        let total_dist: f32 = self.front.iter().map(|n| self.circles[*n].radius).sum();
        match start_index.cmp(&end_index) {
            std::cmp::Ordering::Less => {
                let front_dist: f32 = self.front[start_index..end_index].iter().map(|i| self.circles[*i].radius).sum();
                (
                    2.0 * front_dist,
                    2.0 * total_dist - 2.0 * front_dist - 2.0 * self.circles[self.front[end_index]].radius,
                )
            }
            std::cmp::Ordering::Equal => (0.0, total_dist),
            std::cmp::Ordering::Greater => {
                let back_dist: f32 = self.front[end_index..start_index].iter().map(|i| self.circles[*i].radius).sum();
                (
                    2.0 * total_dist - 2.0 * back_dist - 2.0 * self.circles[self.front[start_index]].radius,
                    2.0 * back_dist - 2.0 * self.circles[self.front[end_index]].radius,
                )
            }
        }
    }
    /// Remove the bits between the start and end and returns the new indexes of the start and end.
    fn front_remove(&mut self, start_index: usize, end_index: usize) -> (usize, usize) {
        if start_index < end_index {
            self.front.drain((start_index + 1)..end_index);
            (start_index, self.front_get_next_index(start_index))
        } else {
            self.front.drain((start_index + 1)..);
            self.front.drain(..end_index);
            (
                self.front.len() - 1,
                0,
            )
        }
    }
}

/// Sorts the radii for you for better results
#[derive(Debug, Clone, Default)]
pub struct CircleSortedPacker {
    radii: Vec<(usize,f32)>,
}

impl CircleSortedPacker {
    pub fn push(&mut self, radius: f32) {
        self.radii.push((self.radii.len(),radius));
        self.radii.sort_by(|a,b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
    }

    pub fn circles(&self) -> Vec<Circle> {
        let mut circles =  vec![Circle::default(); self.radii.len()];
        let mut packer = CirclePacker::default();
        for (_,f) in &self.radii {
            packer.push(*f);
        }
        packer.circles.iter().zip(&self.radii).for_each(|(circle, (i,_))| {
            circles[*i] = circle.clone();
        });
        circles
    }

    pub fn circles_in(&self, embedding_circle: &Circle) -> Vec<Circle> {
        let mut circles =  vec![Circle::default(); self.radii.len()];
        let mut packer = CirclePacker::default();
        for (_,f) in &self.radii {
            packer.push(*f);
        }
        tracing::info!("Packer: {:?}", packer);
        if let Some(radius) = packer.circles.iter().map(|c| c.distance_to_origin() + c.radius).max_by(|a,b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal)) {
            let radius_ratio = embedding_circle.radius / radius;
            tracing::info!("Embedding Radius: {:?}, radius: {}, ratio: {}", embedding_circle.radius, radius, radius_ratio);

            packer.circles.iter().zip(&self.radii).for_each(|(circle, (i,_))| {
                circles[*i].x = circle.x*radius_ratio + embedding_circle.x;
                circles[*i].y = circle.y*radius_ratio + embedding_circle.y;
                circles[*i].radius = circle.radius*radius_ratio;
            });
        }
        tracing::info!("Circles: {:?}", circles);
        circles
    }
}


#[cfg(test)]
mod tests {
    use crate::{Circle, CirclePacker};
    use rand::{prelude::*, rngs::SmallRng};
    use assert_approx_eq::assert_approx_eq;
    fn has_nan(circle: &Circle) -> bool {
        circle.x.is_nan() || circle.y.is_nan() || circle.radius.is_nan()
    }

    fn front_touches(packer: &CirclePacker) {
        tracing::info!("Packer: {:?}", packer);
        for i in 1..packer.front.len() {
            assert_eq!(i,packer.front_get_next_index(i-1));
            let cm = &packer.circles[packer.front[i-1]];
            let cn = &packer.circles[packer.front[i]];
            assert_approx_eq!(cm.distance(cn), cm.radius + cn.radius, 1e-4f32);
        }
        if 1 < packer.front.len() {
            assert_eq!(0,packer.front_get_next_index(packer.front.len() - 1));
            let cm = &packer.circles[packer.front[packer.front.len() - 1]];
            let cn = &packer.circles[packer.front[0]];
            assert_approx_eq!(cm.distance(cn), cm.radius + cn.radius, 1e-4f32);
        }
    }
    #[tracing_test::traced_test]
    #[test]
    pub fn new_node_does_not_intersect_touching() {
        let mut rng = SmallRng::from_seed([0; 32]);
        for _ in 0..20 {
            let mut circle_1 = Circle {
                x: rng.gen(),
                y: rng.gen(),
                radius: rng.gen(),
            };
            
            let mut circle_2 = Circle {
                x: rng.gen(),
                y: rng.gen(),
                radius: 0.0,
            };
            while circle_1.distance(&circle_2) < circle_1.radius {
                circle_1.x = rng.gen();
                circle_1.y = rng.gen();
                circle_1.radius = rng.gen();
                circle_2.x = rng.gen();
                circle_2.y = rng.gen();
            }
            let dist = circle_1.distance(&circle_2);
            circle_2.radius = dist - circle_1.radius;

            let radius: f32 = rng.gen();

            let (circle_3, circle_4) = Circle::new_node(radius, &circle_1, &circle_2);
            assert!(!circle_3.x.is_nan(), "Candidate Circle has a NAN");
            assert!(!circle_3.y.is_nan(), "Candidate Circle has a NAN");
            assert!(!circle_4.x.is_nan(), "Candidate Circle has a NAN");
            assert!(!circle_4.y.is_nan(), "Candidate Circle has a NAN");
            assert_approx_eq!(circle_1.distance(&circle_2), circle_1.radius + circle_2.radius);
            assert_approx_eq!(circle_1.distance(&circle_3), circle_1.radius + circle_3.radius);
            assert_approx_eq!(circle_1.distance(&circle_4), circle_1.radius + circle_4.radius);
            assert_approx_eq!(circle_2.distance(&circle_3), circle_2.radius + circle_3.radius);
            assert_approx_eq!(circle_2.distance(&circle_4), circle_2.radius + circle_4.radius);
            assert!(!circle_1.intersect(&circle_2));
            assert!(!circle_1.intersect(&circle_3));
            assert!(!circle_1.intersect(&circle_4));
            assert!(!circle_2.intersect(&circle_3));
            assert!(!circle_2.intersect(&circle_4));
        }
    }

    fn radii_work(radii: &[f32]) {
        let mut packer = CirclePacker::default();
        for radius in radii {
            let span = tracing::info_span!("Radius: {}", radius);
            let _entered = span.enter();
            packer.push(*radius);
            front_touches(&packer);
            let circles = packer.circles(None);
            for i in 0..circles.len() {
                assert!(!has_nan(&circles[i]), "Circle {} has a nan on radius {}", i, radius);
                for j in 0..i {
                    assert!(!circles[i].intersect(&circles[j]), "{} {:?} intersects with {} {:?} on radius {}", i, &circles[i], j, &circles[j], radius);
                }
                for j in (i+1)..circles.len() {
                    assert!(!circles[i].intersect(&circles[j]), "{} {:?} intersects with {} {:?} on radius {}", i, &circles[i], j, &circles[j], radius);
                }
            }
        }
    }

    #[tracing_test::traced_test]
    #[test]
    pub fn it_works() {
        radii_work(&[50.0, 100.1, 40.2, 30.3, 60.4, 40.5, 80.6, 100.7, 50.8])
    }

    #[tracing_test::traced_test]
    #[test]
    pub fn it_works_2() {
        radii_work(&[32.0, 38.1, 35.2, 1.3, 49.4, 2.5, 85.6])
    }

    #[tracing_test::traced_test]
    #[test]
    pub fn it_works_3() {
        radii_work(&[89.0, 96.1, 8.2, 71.3])
    }

    #[tracing_test::traced_test]
    #[test]
    pub fn it_works_4() {
        radii_work(&[27.0, 56.1, 15.2, 89.3, 91.4, 13.5, 27.6, 86.7])
    }

    #[tracing_test::traced_test]
    #[test]
    pub fn it_works_fuzz() {
        let mut rng = SmallRng::from_seed([0; 32]);
        for _ in 0..20 {
            let mut packer = CirclePacker::default();
            for _ in 0..50 {
                let radius: f32 = rng.gen();
                packer.push(100.0f32 * radius);
                front_touches(&packer);
                let circles = packer.circles(None);
                for i in 0..circles.len() {
                    assert!(!has_nan(&circles[i]), "Circle {} has a nan on radius {}", i, radius);
                    for j in 0..i {
                        assert!(!circles[i].intersect(&circles[j]), "{} {:?} intersects with {} {:?} on radius {}", i, &circles[i], j, &circles[j], radius);
                    }
                    for j in (i+1)..circles.len() {
                        assert!(!circles[i].intersect(&circles[j]), "{} {:?} intersects with {} {:?} on radius {}", i, &circles[i], j, &circles[j], radius);
                    }
                }
            }
        }
    }
}
