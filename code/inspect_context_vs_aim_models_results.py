from imports import *
from util_funcs_2 import *
import model_context_1s
import model_context_ns
import model_context_1f1s
import model_context_1fns
import model_context_nf1s
import model_context_nfns
import model_aim_y


if __name__ == "__main__":
    # NOTE: prep data
    d = load_all_data_grand()

    d = d[d.group.isin([7, 8, 19])]

    d["group_class"] = d["group"].map({7: "blocked", 8: "blocked", 19: "interleaved"})

    d["group"] = d["group"].astype("category")
    d["sig_mpep_prev"] = d["sig_mpep_prev"].astype("category")

    d["group"] = d["group"].cat.remove_unused_categories()
    d["phase"] = d["phase"].cat.remove_unused_categories()
    d["sig_mpep_prev"] = d["sig_mpep_prev"].cat.remove_unused_categories()

    # dd = d[d.phase.isin(["adaptation", "washout"])]
    dd = d[d.phase.isin(["adaptation"])]
    # dd = dd[dd.group.isin([19])]

    froot = "../fits/adaptation_only"

    # models = ["1s", "ns", "1f1s", "1fns", "nf1s", "nfns", "aim_y"]
    models = ["1s", "ns", "aim_y"]

    d_subject_list = []
    d_params_list = []
    for file in os.listdir(froot):
        if file.endswith(".txt"):
            subject = file.split("_")[-1].split(".")[0]
            group = file.split("_")[-3]
            subject = int(subject)
            group = int(group)

            d_subject = dd[(dd["subject"] == subject) & (dd["group"] == group)].copy()

            if not d_subject.empty:
                trial = d_subject.trial.to_numpy()
                rot = d_subject.rot.to_numpy()
                sig_mp = d_subject.sig_mp.to_numpy()
                sig_ep = d_subject.sig_ep.to_numpy()
                sig_mpep_prev = d_subject.sig_mpep_prev.to_numpy()
                group = d_subject.group.to_numpy()
                x_obs_mp = d_subject["ha_init"].to_numpy()
                x_obs_ep = d_subject["ha_end"].to_numpy()

                args = (rot, sig_mp, sig_ep, x_obs_mp, x_obs_ep)
                params = np.loadtxt(os.path.join(froot, file), delimiter=",")[:-2]

                p_names = np.array(
                    [
                        "alpha_ff",
                        "beta_ff",
                        "bias_ff",
                        "alpha_ff2",
                        "beta_ff2",
                        "bias_ff2",
                        "alpha_fb",
                        "beta_fb",
                        "xfb_init",
                        "gamma_fbint_1",
                        "gamma_fbint_2",
                        "gamma_fbint_3",
                        "gamma_fbint_4",
                        "gamma_ff_1",
                        "gamma_ff_2",
                        "gamma_ff_3",
                        "gamma_ff_4",
                    ]
                )

                d_params = pd.DataFrame({"params": params, "p_names": p_names})
                d_params["subject"] = np.unique(subject)[0]
                d_params["group"] = np.unique(group)[0]

                mod_name = ""
                if "context_1s" in file:
                    mod_name = "1s"
                    k = 18 - 3
                    x_pred = model_context_1s.simulate(params, args)

                if "context_ns" in file:
                    mod_name = "ns"
                    k = 18 - 3
                    x_pred = model_context_ns.simulate(params, args)

                if "context_1f1s" in file:
                    mod_name = "1f1s"
                    k = 18
                    x_pred = model_context_1f1s.simulate(params, args)

                if "context_1fns" in file:
                    mod_name = "1fns"
                    k = 18
                    x_pred = model_context_1fns.simulate(params, args)

                if "context_nf1s" in file:
                    mod_name = "nf1s"
                    k = 18
                    x_pred = model_context_nf1s.simulate(params, args)

                if "context_nfns" in file:
                    mod_name = "nfns"
                    k = 18
                    x_pred = model_context_nfns.simulate(params, args)

                if "aim_y" in file:
                    mod_name = "aim_y"
                    k = 18 - 3
                    x_pred = model_aim_y.simulate(params, args)

                if mod_name == "":
                    print(file)

                x_pred_mp = x_pred[1]
                x_pred_ep = x_pred[0]

                ss_tot_mp = np.nansum((x_obs_mp - np.nanmean(x_obs_mp)) ** 2)
                ss_reg_mp = np.nansum((x_pred_mp - np.nanmean(x_pred_mp)) ** 2)
                ss_res_mp = np.nansum((x_obs_mp - x_pred_mp) ** 2)

                ss_tot_ep = np.nansum((x_obs_ep - np.nanmean(x_obs_ep)) ** 2)
                ss_reg_ep = np.nansum((x_pred_ep - np.nanmean(x_pred_ep)) ** 2)
                ss_res_ep = np.nansum((x_obs_ep - x_pred_ep) ** 2)

                r_squared_mp = 1 - ss_res_mp / ss_tot_mp
                r_squared_ep = 1 - ss_res_ep / ss_tot_ep
                r_squared = 1 - (ss_res_ep + ss_res_mp) / (ss_tot_ep + ss_tot_mp)

                n = d_subject.shape[0]
                bic = n * np.log(1 - r_squared) + k * np.log(n)

                d_subject["r_squared"] = r_squared
                d_subject["bic"] = bic
                d_subject["ha_init_pred"] = x_pred_mp
                d_subject["ha_end_pred"] = x_pred_ep
                d_subject["model"] = mod_name

                d_params["model"] = mod_name

                d_subject_list.append(d_subject)
                d_params_list.append(d_params)

    d_subject = pd.concat(d_subject_list)
    d_params = pd.concat(d_params_list)

    # NOTE: plot fits and  best fitting parameters
    bounds = (
        # alpha_ff, beta_ff, bias_ff,
        (0, 1),
        (0, 1),
        (-20, 20),
        # alpha_ff2, beta_ff2, bias_ff2,
        (0, 1),
        (0, 1),
        (0, 0),  # NOTE: bias terms are additive
        # alpha_fb, beta_fb, xfb_init,
        (0, 1),
        (0, 1),  # NOTE: This was (-10, 10)
        (-2, 2),
        # gamma_fbint_1, gamma_fbint_2, gamma_fbint_3, gamma_fbint_4,
        (0, 1),
        (0, 1),
        (0, 1),
        (0, 1),
        # gamma_ff_1, gamma_ff_2, gamma_ff_3, gamma_ff_4,
        (0, 1),
        (0, 1),
        (0, 1),
        (0, 1),
        # temporal_discount,
        (0, 0),  # NOTE: td doesn't do anything with matched mp/ep
    )

    for m in models:
        bounds_2 = list(bounds)
        if m == "1s":
            bounds_2[3] = (1, 1)
            bounds_2[4] = (0, 0)
            bounds_2[5] = (0, 0)
        if m == "ns":
            bounds_2[3] = (1, 1)
            bounds_2[4] = (0, 0)
            bounds_2[5] = (0, 0)
        if m == "aim_y":
            bounds_2[3] = (1, 1)
            bounds_2[4] = (0, 0)
            bounds_2[5] = (0, 0)
            bounds_2[-2] = (-20, 20)
            bounds_2[-3] = (-20, 20)
            bounds_2[-4] = (-20, 20)
            bounds_2[-5] = (-20, 20)

        fig = plt.figure(figsize=(14, 5))
        gs = fig.add_gridspec(2, 17)
        gs.update(wspace=2.5, bottom=0.2)

        ax = fig.add_subplot(gs[0, :])
        sns.lineplot(
            # data=d_subject[(d_subject.model == m)],
            data=d_subject[(d_subject.model == m) & (d_subject.phase == "adaptation")],
            x="trial",
            y="ha_init_pred",
            hue="sig_mpep_prev",
            style="sig_mpep_prev",
            markers=True,
            legend=False,
            ax=ax,
        )
        ax.set_ylim([0, 12])

        for i, p in enumerate(d_params.p_names.unique()):
            ax = fig.add_subplot(gs[1, i])
            sns.boxplot(
                data=d_params[(d_params.model == m) & (d_params.p_names == p)],
                x="p_names",
                y="params",
                ax=ax,
                whis=(0, 100),
            )
            sns.scatterplot(
                data=d_params[(d_params.model == m) & (d_params.p_names == p)],
                x="p_names",
                y="params",
                ax=ax,
            )
            ax.set_ylim(bounds_2[i])
            ax.set_yticks(bounds_2[i])
            ax.set_ylabel("")
            ax.set_xlabel("")
            for tick in ax.get_xticklabels():
                tick.set_rotation(45)

        plt.suptitle(m)
        plt.savefig(f"../figures/fig_{m}_params.pdf")
        plt.close()

    # NOTE: Plot fits vs observed
    fig, ax = plt.subplots(3, 2, squeeze=False, figsize=(8, 11))

    ds = d_subject[
        ["subject", "phase", "trial", "ha_init", "sig_mpep_prev"]
    ].drop_duplicates()
    sns.lineplot(
        data=ds,
        x="trial",
        y="ha_init",
        hue="sig_mpep_prev",
        style="phase",
        markers=True,
        legend=False,
        ax=ax[0, 0],
    )
    ax[0, 0].set_title("Human Data")

    ds = d_subject[d_subject["model"] == "1s"]
    ds = ds[ds["trial"] < ds["trial"].max()]
    sns.lineplot(
        data=ds,
        x="trial",
        y="ha_init_pred",
        hue="sig_mpep_prev",
        style="phase",
        markers=True,
        legend=False,
        ax=ax[0, 1],
    )
    ax[0, 1].set_title("Standard Model")

    ds = d_subject[d_subject["model"] == "ns"]
    ds = ds[ds["trial"] < ds["trial"].max()]
    sns.lineplot(
        data=ds,
        x="trial",
        y="ha_init_pred",
        hue="sig_mpep_prev",
        style="phase",
        markers=True,
        legend=False,
        ax=ax[1, 0],
    )
    ax[1, 0].set_title("Context Model")

    ds = d_subject[d_subject["model"] == "aim_y"]
    ds = ds[ds["trial"] < ds["trial"].max()]
    sns.lineplot(
        data=ds,
        x="trial",
        y="ha_init_pred",
        hue="sig_mpep_prev",
        style="phase",
        markers=True,
        legend=False,
        ax=ax[1, 1],
    )
    ax[1, 1].set_title("Aiming Model")

    d = d_subject[["group", "subject", "model", "bic"]].drop_duplicates().copy()
    d["model"] = d["model"].replace(
        {"1s": "Standard", "ns": "Context", "aim_y": "Aiming"}
    )
    d["model"] = d["model"].astype("category")
    d["model"] = d["model"].cat.reorder_categories(["Standard", "Context", "Aiming"])
    sns.boxplot(data=d, x="model", y="bic", color=".8", whis=(0, 100), ax=ax[2, 0])
    sns.stripplot(data=d, x="model", y="bic", color=".3", ax=ax[2, 0])
    ax[2, 0].set_title("Model BICs")

    labels = ["a", "b", "c", "d", "e", ""]
    for i, axx in enumerate(ax.flatten()):
        label = labels[i]
        axx.text(-0.1, 1.1, label, transform=axx.transAxes, fontsize=16)

    ax[-1, -1].remove()

    [x.set_ylabel("Hand Angle (deg)") for x in ax.flatten()[:-2]]
    [x.set_ylim([0, 14]) for x in ax.flatten()[:-2]]
    plt.tight_layout()
    plt.savefig("../figures/fig_obs_vs_pred.pdf")
    plt.savefig("../figures/fig1.tif")
    plt.savefig("../figures/fig1.png")
    plt.close()

    # NOTE: model selection
    # d = d_subject[["group", "subject", "model", "bic"]].drop_duplicates().copy()
    # d["model"] = d["model"].replace({"1s": "Standard", "ns": "Context", "aim_y": "Aiming"})
    # d["model"] = d["model"].astype("category")
    # d["model"] = d["model"].cat.reorder_categories(["Standard", "Context", "Aiming"])

    # fig, ax = plt.subplots(1, 1, squeeze=False, figsize=(6, 6))
    # ax[0, 0].set_title("Model BICs")
    # sns.boxplot(data=d, x="model", y="bic", color=".8", ax=ax[0, 0], whis=(0, 100))
    # sns.stripplot(data=d, x="model", y="bic", ax=ax[0, 0], color=".3")
    # plt.savefig("../figures/fig_model_bic.pdf")
    # plt.savefig("../figures/fig2.tif")
    # plt.close()

    # def rank_models(x):
    #     x = x.sort_values(by="bic")
    #     x["rank"] = np.arange(1, x.shape[0] + 1, 1)
    #     return x

    # rank_table = []
    # for g in [19]:
    #     dd = d[d["group"] == g]
    #     dd = dd.groupby(["subject"]).apply(rank_models).reset_index(drop=True)
    #     rank_table.append(pd.crosstab(dd["model"], dd["rank"]))

    # fig, ax = plt.subplots(1, 1, squeeze=False, figsize=(6, 6))
    # cm = sns.color_palette("Blues", as_cmap=True)
    # sns.heatmap(
    #     rank_table[0],
    #     ax=ax[0, 0],
    #     annot=True,
    #     cbar=False,
    #     cmap=sns.color_palette("Blues", as_cmap=True),
    # )
    # ax[0, 0].set_title("Model Ranks")
    # plt.tight_layout()
    # plt.savefig("../figures/fig_model_rank.pdf")
    # plt.close()

    # NOTE: figure for ACNS 2023 poster
    # NOTE: figure for ACNS 2023 poster
    # NOTE: figure for ACNS 2023 poster
    # NOTE: figure for ACNS 2023 poster
    fig, ax = plt.subplots(1, 3, squeeze=False, figsize=(12, 4))

    ds = d_subject[d_subject["model"] == "1s"]
    ds = ds[ds["trial"] < ds["trial"].max()]
    sns.lineplot(
        data=ds,
        x="trial",
        y="ha_init_pred",
        hue="sig_mpep_prev",
        style="phase",
        markers=True,
        legend=False,
        ax=ax[0, 0],
    )
    ax[0, 0].set_title("Standard Model")

    ds = d_subject[d_subject["model"] == "ns"]
    ds = ds[ds["trial"] < ds["trial"].max()]
    sns.lineplot(
        data=ds,
        x="trial",
        y="ha_init_pred",
        hue="sig_mpep_prev",
        style="phase",
        markers=True,
        legend=False,
        ax=ax[0, 1],
    )
    ax[0, 1].set_title("Context Model")

    ds = d_subject[d_subject["model"] == "aim_y"]
    ds = ds[ds["trial"] < ds["trial"].max()]
    sns.lineplot(
        data=ds,
        x="trial",
        y="ha_init_pred",
        hue="sig_mpep_prev",
        style="phase",
        markers=True,
        legend=False,
        ax=ax[0, 2],
    )
    ax[0, 2].set_title("Aiming Model")

    labels = ["a", "b", "c"]
    for i, axx in enumerate(ax.flatten()):
        label = labels[i]
        axx.text(-0.1, 1.1, label, transform=axx.transAxes, fontsize=16)

    [x.set_ylabel("Hand Angle (deg)") for x in ax.flatten()[:-1]]
    [x.set_ylim([0, 14]) for x in ax.flatten()[:-1]]
    plt.tight_layout()
    plt.savefig("../figures/fig1_context_vs_aim_poster_1.png")
    plt.close()


    # d = d_subject[["group", "subject", "model", "bic"]].drop_duplicates().copy()
    # d["model"] = d["model"].replace(
    #     {"1s": "Standard", "ns": "Context", "aim_y": "Aiming"}
    # )
    # d["model"] = d["model"].astype("category")
    # d["model"] = d["model"].cat.reorder_categories(["Standard", "Context", "Aiming"])
    # sns.boxplot(data=d, x="model", y="bic", color=".8", whis=(0, 100), ax=ax[1, 1])
    # sns.stripplot(data=d, x="model", y="bic", color=".3", ax=ax[1, 1])
    # ax[1, 1].set_title("Model BICs")
