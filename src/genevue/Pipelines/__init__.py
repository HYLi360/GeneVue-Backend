#  Copyright (C) 2025-2026, HYLi360.
#  Free software distributed under the terms of the GNU GPL-3.0 license,
#  and comes with ABSOLUTELY NO WARRANTY.
#  See at <https://www.gnu.org/licenses/gpl-3.0.en.html>

"""
A series of pre-defined pipeline.
"""

import inspect
from collections import defaultdict, deque
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Callable, Dict, List, Optional, get_type_hints


class NodeStatus(Enum):
    PENDING = "pending"
    RUNNING = "running"
    SUCCESS = "success"
    FAILED = "failed"
    SKIPPED = "skipped"


@dataclass
class ProcessNode:
    func: Callable
    outputs: List[str] = field(default_factory=list)
    depends: Optional[List[str]] = None
    status: NodeStatus = NodeStatus.PENDING
    results: Dict[str, Any] = field(default_factory=dict, repr=False)


class Pipeline:
    def __init__(self):
        self.nodes: Dict[str, ProcessNode] = {}
        self.inputs: Dict[str, Any] = {}

    def node(self, node_name: str, outputs=None, depends=None):
        if outputs is None:
            outputs = []

        def decorator(func):
            self.nodes[node_name] = ProcessNode(
                func=func,
                outputs=outputs,
                depends=depends or [],
            )
            return func

        return decorator

    def _build(self):
        indegree = defaultdict(int)
        adj_list = defaultdict(list)

        for node_name, node in self.nodes.items():
            for depend in node.depends:
                if depend not in self.nodes:
                    raise ValueError(
                        f"Node {node_name} depends on unknown node {depend}"
                    )
                indegree[node_name] += 1
                adj_list[depend].append(node_name)

        execute_order: List[str] = []
        queue = deque([node for node in self.nodes if indegree[node] == 0])

        while queue:
            current = queue.popleft()
            execute_order.append(current)

            for next_node in adj_list[current]:
                indegree[next_node] -= 1
                if indegree[next_node] == 0:
                    queue.append(next_node)

        if len(execute_order) != len(self.nodes):
            raise ValueError("Loop dependency detected!")

        return execute_order

    def _collect_inputs(self, node: ProcessNode) -> Dict[str, Any]:
        sig = inspect.signature(node.func)
        args_pool = self.inputs
        kwargs = {}

        if node.depends:
            for depend in node.depends:
                depended_node = self.nodes[depend]
                for result in depended_node.results:
                    args_pool[result] = depended_node.results[result]

        for name, param in sig.parameters.items():
            if name not in args_pool:
                raise KeyError(f"Missing input: {name}")
            kwargs[name] = self.inputs[name]

        return kwargs

    def _bind_outputs(self, node_name: str, outputs: List[str], result: Any):
        if not outputs:
            return

        if isinstance(result, dict):
            for name in outputs:
                if name not in result:
                    raise KeyError(f"Missing output key: {name}")
                self.nodes[node_name].results[name] = result[name]
            return

        if not isinstance(result, (tuple, list)):
            raise TypeError(
                f"Node has multiple outputs {outputs}, but return value is not tuple/list"
            )

        if len(outputs) == 1:
            self.nodes[node_name].results[outputs[0]] = result
            return

        if len(result) != len(outputs):
            print(result, outputs)
            raise ValueError(
                f"Output count mismatch: expected {len(outputs)}, got {len(result)}"
            )

        for name, value in zip(outputs, result):
            self.nodes[node_name].results[name] = value

    def check_needed_inputs(self):
        produced, consumed = [], []
        for node in self.nodes.values():
            produced.extend(node.outputs)
            consumed.extend(get_type_hints(node.func).keys())

        res = list(set(consumed) - set(produced))
        res.sort()

        print(res)

    def execute(self, input_data=None):
        self.inputs = dict(input_data or {})
        execute_order = self._build()

        for node_name in execute_order:
            node = self.nodes[node_name]
            node.status = NodeStatus.RUNNING
            kwargs = self._collect_inputs(node)
            result = node.func(**kwargs)
            self._bind_outputs(node_name, node.outputs, result)
            node.status = NodeStatus.SUCCESS

    def dry_run(self):
        execute_order = self._build()
        print(execute_order)
