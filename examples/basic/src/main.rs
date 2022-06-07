use yew::prelude::*;
use circle_packer::{CirclePacker, Circle};

const DEPTH_COLOR: [&'static str;4] = ["blue", "red", "green", "yellow"];

struct CircleNode;

impl Component for CircleNode {
    type Message = ();
    type Properties = CircleNodeProps;

    fn create(_ctx: &Context<Self>) -> Self {
        CircleNode
    }
    fn view(&self, ctx: &Context<Self>) -> Html {
        let container = ctx.props().contain.clone();
        let mut packer = CirclePacker::default();
        ctx.props().nodes.iter().for_each(|n| packer.push(n.radius));
        let circles = packer.circles_in(&container);
        let color = DEPTH_COLOR[ctx.props().depth];
        
        html! {
            <g>
                <circle cx={container.x.to_string()} cy={container.y.to_string()} r={container.radius.to_string()} fill={color}/>
                {
                    ctx.props().nodes.iter().zip(circles).map(|(node, circle)|  {
                        html! {
                            <CircleNode nodes={node.descendants.clone()} contain={circle} depth={ctx.props().depth + 1}/>
                        }
                    }).collect::<Html>()
                }
            </g>
        }
    }
}

#[derive(PartialEq, Properties, Default)]
struct CircleNodeProps {
    depth: usize,
    contain: Circle,
    nodes: Vec<Node>,
}

#[derive(PartialEq, Default, Clone)]
struct Node {
    radius: f32,
    descendants: Vec<Node>,
}

struct CircleTree;

#[derive(PartialEq, Properties, Default)]
struct CircleTreeProps {
    root: Node,
} 

impl Component for CircleTree {
    type Message = ();
    type Properties = CircleTreeProps;

    fn create(_ctx: &Context<Self>) -> Self {
        CircleTree
    }
    fn view(&self, ctx: &Context<Self>) -> Html {
        let root = &ctx.props().root;
        let root_circle = Circle {
            x: 0.0,
            y: 0.0,
            radius: root.radius,
        };
        
        html! {
            <svg viewBox="-500 -500 1000 1000" xmlns="http://www.w3.org/2000/svg">
                <CircleNode nodes={root.descendants.clone()} contain={root_circle} depth={0}/>
            </svg>
        }
    }
}

pub fn main() {
    let root = Node {radius: 350.0, descendants: vec![
        Node {radius: 38.1, descendants: vec![]},
        Node {radius: 35.2, descendants: vec![]},
        Node {radius: 1.3, descendants: vec![]},
        Node {radius: 49.4, descendants: vec![]},
        Node {radius: 2.5, descendants: vec![]},
        Node {radius: 85.6, descendants: vec![
            Node {radius: 32.0, descendants: vec![]},
            Node {radius: 38.1, descendants: vec![]},
            Node {radius: 35.2, descendants: vec![]},
            Node {radius: 1.3, descendants: vec![]},
            Node {radius: 49.4, descendants: vec![]},
            Node {radius: 2.5, descendants: vec![]},
            Node {radius: 85.6, descendants: vec![
                Node {radius: 38.1, descendants: vec![]},
                Node {radius: 35.2, descendants: vec![]},
                Node {radius: 1.3, descendants: vec![]},
                Node {radius: 49.4, descendants: vec![]},
            ]},
        ]},
        Node {radius: 84.7, descendants: vec![
            Node {radius: 32.0, descendants: vec![]},
            Node {radius: 38.1, descendants: vec![]},
            Node {radius: 35.2, descendants: vec![]},
        ]},
        Node {radius: 29.8, descendants: vec![]},
        Node {radius: 7.9, descendants: vec![]},
        Node {radius: 31.10, descendants: vec![]},
        Node {radius: 6.11, descendants: vec![]},
        Node {radius: 11.12, descendants: vec![]},
    ]};
    tracing_wasm::set_as_global_default();
    yew::start_app_with_props::<CircleTree>(CircleTreeProps {root});
}
