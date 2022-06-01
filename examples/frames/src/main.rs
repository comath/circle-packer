use yew::prelude::*;
use yew_circle_packing::{CirclePacker, Circle};

const DEPTH_COLOR: [&'static str;3] = ["blue", "red", ];

#[derive(PartialEq, Properties, Default)]
pub struct CircleProperies {
    circle: Circle,
    index: usize,
    color: String,
}

pub struct CircleComponent;

impl Component for CircleComponent {
    type Message = ();
    type Properties = CircleProperies;

    fn create(_ctx: &Context<Self>) -> Self {
        CircleComponent
    }
    fn view(&self, ctx: &Context<Self>) -> Html {
        html! {
            <g>
                <circle cx={ctx.props().circle.x.to_string()} cy={ctx.props().circle.y.to_string()} r={ctx.props().circle.radius.to_string()} fill-color={ctx.props().color}/>
                <g>
                </g>
            </g>
        }
    }
}

pub struct CircleDemoFrame;

#[derive(PartialEq, Properties)]
pub struct CircleDemoFrameProps {
    circles: Vec<Circle>,
}

impl Component for CircleDemoFrame {
    type Message = ();
    type Properties = CircleDemoFrameProps;
    fn create(_ctx: &Context<Self>) -> Self {
        CircleDemoFrame
    }
    
    fn view(&self, ctx: &Context<Self>) -> Html {
        let circles = &ctx.props().circles;
        html! {
            <svg viewBox="-500 -500 1000 1000" xmlns="http://www.w3.org/2000/svg">
                {
                    circles.into_iter().enumerate().map(|(index,circle)| {
                        html! {
                            <CircleComponent circle={circle.clone()} {index} opacity={1.0}/>
                        }
                    }).collect::<Html>()
                }
                <circle cx="0" cy="0" r="4" fill="red" />
            </svg>
        }
    }
}

pub struct CircleDemo;

#[derive(PartialEq, Properties)]
pub struct CircleDemoProps {
    radii: Vec<f32>,
}

impl Component for CircleDemo {
    type Message = ();
    type Properties = CircleDemoProps;
    fn create(_ctx: &Context<Self>) -> Self {
        CircleDemo
    }
    fn view(&self, ctx: &Context<Self>) -> Html {
        let mut packer = CirclePacker::default();
        packer.push(ctx.props().radii[0]);
        let radii = ctx.props().radii[1..].to_vec();

        html! {
            <div>
                {
                    radii.into_iter().map(|radius| {
                        
                        let circles= packer.circles();
                        packer.push(radius);

                        html! {
                            <CircleDemoFrame {circles} />
                        }
                    }).collect::<Html>()
                }
                <CircleDemoFrame circles={packer.circles()} />
            </div>
        }
    }
}

pub fn main() {
    tracing_wasm::set_as_global_default();
    yew::start_app_with_props::<CircleDemo>(CircleDemoProps {
        radii: vec![32.0, 38.1, 35.2, 1.3, 49.4, 2.5, 85.6, 84.7, 29.8, 7.9, 31.10, 6.11, 11.12],
    });
}
