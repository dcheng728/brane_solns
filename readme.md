# On the 12d Interpretation of type IIB String Theory
- The primary (scientific) objective is to explore whether there exists a 12d interpretation of the type IIB string theory. 
- The secondary objective is to explore how AI can be helpful in theoretical physics research.

# How will AI be used?
Modern ML has been applied across many branches of science, especially those dealing with numerical data, since such data is naturally amenable to existing optimization and ML techniques. Theoretical physics and pure math, by contrast, revolve around abstract ideas and symbolic computation, so the role of AI in these fields has historically been more limited.

Language models may change this. They appear capable of manipulating abstract concepts and performing symbolic computations in ways that go beyond mechanical pattern-matching, which makes the question of AI's role in theoretical physics worth revisiting. Other physicists have already found LLMs useful for tasks like literature summaries and sanity-checking understanding through browser-based conversations — but I believe their utility need not be restricted to that of a chatbot.

This repo exists to explore that possibility. Alongside my master's dissertation, I want to investigate how AI can facilitate theoretical physics research — both the creative aspects of the work and the more menial tasks that AI could take off a physicist's hands.

To keep the scope manageable and the results trustworthy, an agentic system seems like the right approach. For instance, rather than having a single AI try to solve a given equation of motion, one can set up two agents - a proposer and a verifier - in an adversarial loop, where the proposer suggests a solution and the verifier checks it both symbolically and numerically. The proposer might be granted read-only access, while the verifier holds both read and write privileges.

# Organization


to start:

```
conda env create -f environment.yml
```