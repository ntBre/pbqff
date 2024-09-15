async function qffbuddy(self) {
	console.log(self.geometry.value)
	let optimize = document.querySelector("[name=optimize]").checked;
	let findiff = document.querySelector("[name=findiff]").checked;
	console.log(optimize);

	let queue_template;
	if (self.queue_template.value) {
		queue_template = `queue_template = """
${self.queue_template.value}
"""`;
	} else {
		queue_template = "";
	}
	let hybrid_template;
	if (self.hybrid_template.value) {
		hybrid_template = `hybrid_template = """
${self.hybrid_template.value}
"""`;
	} else {
		hybrid_template = "";
	}
	let template = `\
geometry = """
${self.geometry.value}
"""
optimize = ${optimize}
charge = ${self.charge.value}
step_size = ${self.step_size.value}
sleep_int = ${self.sleep_int.value}
job_limit = ${self.job_limit.value}
chunk_size = ${self.chunk_size.value}
check_int = ${self.check_int.value}
coord_type = "${self.coord_type.value}"
findiff = ${findiff}
program = "${self.program.value}"
queue = "${self.queue.value}"
template = """
${self.template.value}
"""
${queue_template}
${hybrid_template}
`;
	try {
		await navigator.clipboard.writeText(template);
	} catch (error) {
		console.error(error.message);
	}
	return false;
}
